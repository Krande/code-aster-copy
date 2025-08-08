! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine nmdopi(modelz, numedd, ds_algopara, sdpilo)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/nueqch.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
!
    character(len=*), intent(in) :: modelz
    character(len=24), intent(in) :: numedd
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=19), intent(in) :: sdpilo
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (STRUCTURES DE DONNEES)
!
! CONSTRUCTION DE LA SD PILOTAGE
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL
! In  ds_algopara      : datastructure for algorithm parameters
! OUT SDPILO : SD PILOTAGE
!               .PLTK
!                (1) = TYPE DE PILOTAGE
!                (2) = LIGREL POUR LES PILOTAGES PAR ELEMENTS
!                (3) = NOM DE LA CARTE DU TYPE (PILO_K)
!                (4) = NOM DE LA CARTE DU TYPE (PILO_R) MIN/MAX
!                (5) = PROJECTION 'OUI' OU 'NON' SUR LES BORNES
!                (6) = TYPE DE SELECTION : 'RESIDU',
!                        'NORM_INCR_DEPL' OU 'ANGL_INCR_DEPL'
!                (7) = EVOLUTION DES BORNES
!                        'CROISSANT', 'DECROISSANT' OU 'SANS'
!               .PLCR  COEFFICIENTS DU PILOTAGE
!               .PLIR  PARAMETRES DU PILOTAGE
!                (1) = COEF_PILO
!                (2) = ETA_PILO_MAX
!                (3) = ETA_PILO_MIN
!                (4) = ETA_PILO_R_MAX
!                (5) = ETA_PILO_R_MIN
!                (6) = COEF_PILO AU PAS DE TEMPS CONVERGE PRECEDENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbno, numequ, nddl, nb_node_mesh, nb_dof_acti
    integer(kind=8) :: nume_node
    integer(kind=8) :: ino, iddl
    integer(kind=8) :: jvale
    integer(kind=8) :: jplir, jpltk
    integer(kind=8) :: ibid, n1, n2, neq, ndim
    real(kind=8) :: coef, lm(2)
    character(len=8) :: mesh, lborn(2), nomcmp
    character(len=8) :: modele
    character(len=16) :: relmet
    character(len=24) :: lisnoe, liscmp
    integer(kind=8) :: jlinoe, jlicmp
    character(len=24) :: typpil, projbo, typsel, evolpa, txt(2)
    character(len=19) :: chapil, selpil, ligrmo, ligrpi
    character(len=19) :: careta, cartyp
    real(kind=8) :: etrmax, etrmin, etamin, etamax
    integer(kind=8) :: nbmocl
    character(len=16) :: limocl(2), tymocl(2)
    integer(kind=8) :: ifm, niv
    real(kind=8), pointer :: plsl(:) => null()
    aster_logical :: lSelectDof
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_17')
    end if

! - INITIALISATIONS
    modele = modelz
    call dismoi('NOM_MAILLA', numedd, 'NUME_DDL', repk=mesh)
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_node_mesh)
    call dismoi('DIM_GEOM', mesh, 'MAILLAGE', repi=ndim)
    lisnoe = '&&NMDOPI.LISNOE'
    liscmp = '&&NMDOPI.LISCMP'
    nbmocl = 2
    limocl(1) = 'GROUP_NO'
    limocl(2) = 'NOEUD'
    tymocl(1) = 'GROUP_NO'
    tymocl(2) = 'NOEUD'
    nbno = 0

! --- LECTURE DU TYPE ET DE LA ZONE
    call wkvect(sdpilo(1:19)//'.PLTK', 'V V K24', 7, jpltk)
    call getvtx('PILOTAGE', 'TYPE', iocc=1, scal=typpil, nbret=n1)
    zk24(jpltk) = typpil
    call getvtx('PILOTAGE', 'PROJ_BORNES', iocc=1, scal=projbo, nbret=n1)
    zk24(jpltk+4) = projbo
    call getvtx('PILOTAGE', 'SELECTION', iocc=1, scal=typsel, nbret=n1)
    zk24(jpltk+5) = typsel
    call getvtx('PILOTAGE', 'EVOL_PARA', iocc=1, scal=evolpa, nbret=n1)
    zk24(jpltk+6) = evolpa

! - PARAMETRES COEF_MULT ET ETA_PILO_MAX
    call wkvect(sdpilo(1:19)//'.PLIR', 'V V R8', 6, jplir)
    call getvr8('PILOTAGE', 'COEF_MULT', iocc=1, scal=coef, nbret=n1)
    zr(jplir) = coef
    zr(jplir+5) = coef
    if (abs(coef) .le. r8prem()) then
        call utmess('F', 'PILOTAGE_3')
    end if
!
    call getvr8('PILOTAGE', 'ETA_PILO_R_MAX', iocc=1, scal=etrmax, nbret=n1)
    if (n1 .ne. 1) etrmax = r8gaem()
    zr(jplir+3) = etrmax
!
    call getvr8('PILOTAGE', 'ETA_PILO_R_MIN', iocc=1, scal=etrmin, nbret=n2)
    if (n2 .ne. 1) etrmin = -r8gaem()
    zr(jplir+4) = etrmin
!
    call getvr8('PILOTAGE', 'ETA_PILO_MAX', iocc=1, scal=etamax, nbret=n1)
    if (n1 .ne. 1) then
        etamax = r8vide()
    else
        if (etamax .gt. zr(jplir+3)) then
            call utmess('F', 'PILOTAGE_48')
        end if
    end if
    zr(jplir+1) = etamax
!
    call getvr8('PILOTAGE', 'ETA_PILO_MIN', iocc=1, scal=etamin, nbret=n2)
    if (n2 .ne. 1) then
        etamin = r8vide()
    else
        if (etamin .lt. zr(jplir+4)) then
            call utmess('F', 'PILOTAGE_49')
        end if
    end if
    zr(jplir+2) = etamin
!

! ======================================================================
!             PILOTAGE PAR PREDICTION ELASTIQUE : PRED_ELAS
! ======================================================================
!
    if (typpil .eq. 'PRED_ELAS' .or. typpil .eq. 'DEFORMATION') then
!
        call exlima('PILOTAGE', 1, 'V', modele, ligrpi)
        zk24(jpltk+1) = ligrpi
!
!
        cartyp = '&&NMDOPI.TYPEPILO'
        call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrmo)
        call mecact('V', cartyp, 'MODELE', ligrmo, 'PILO_K', &
                    ncmp=1, nomcmp='TYPE', sk=typpil)
        zk24(jpltk+2) = cartyp
!
        lm(1) = etrmax
        lm(2) = etrmin
        careta = '&&NMDOPI.BORNEPILO'
        lborn(1) = 'A0'
        lborn(2) = 'A1'
        call mecact('V', careta, 'MODELE', ligrmo, 'PILO_R', &
                    ncmp=2, lnomcmp=lborn, vr=lm)
        zk24(jpltk+3) = careta
!
!
!
! ======================================================================
!              PILOTAGE PAR UN DEGRE DE LIBERTE : DDL_IMPO
! ======================================================================
!
    else if (typpil .eq. 'DDL_IMPO') then
!
        call reliem(modele, mesh, 'NU_NOEUD', 'PILOTAGE', 1, &
                    nbmocl, limocl, tymocl, lisnoe, nbno)
        if (typpil .eq. 'DDL_IMPO') then
            if (nbno .ne. 1) then
                call utmess('F', 'PILOTAGE_50')
            end if
            coef = 1.d0
        end if
!
!
! ======================================================================
!      PILOTAGE PAR UNE METHODE DE TYPE LONGUEUR D'ARC : LONG_ARC
! ======================================================================
!
    else if (typpil .eq. 'LONG_ARC') then
!
        call reliem(modele, mesh, 'NU_NOEUD', 'PILOTAGE', 1, &
                    nbmocl, limocl, tymocl, lisnoe, nbno)
        if (typpil .eq. 'LONG_ARC') then
            if (nbno .eq. 0) then
                call utmess('F', 'PILOTAGE_57')
            end if
            coef = 1.d0/nbno
        end if
    end if
!
! --- CREATION SD SELECTION DES DDLS EN FEM ?
!
    lSelectDof = ((typpil .eq. 'LONG_ARC') .or. (typpil .eq. 'DDL_IMPO'))

    if (lSelectDof) then
        call getvtx('PILOTAGE', 'NOM_CMP', iocc=1, nbval=0, nbret=nddl)
        nddl = -nddl
        if (nddl .ne. 1 .and. typpil .eq. 'DDL_IMPO') then
            txt(1) = 'NOM_CMP'
            txt(2) = typpil
            call utmess('F', 'PILOTAGE_56', nk=2, valk=txt)
        else if (nddl .eq. 0 .and. typpil .eq. 'LONG_ARC') then
            txt(1) = 'NOM_CMP'
            txt(2) = typpil
            call utmess('F', 'PILOTAGE_55', nk=2, valk=txt)
        end if
        if (nddl .gt. 0) then
            call wkvect(liscmp, 'V V K8', nddl, jlicmp)
            call getvtx('PILOTAGE', 'NOM_CMP', iocc=1, nbval=nddl, vect=zk8(jlicmp), &
                        nbret=ibid)
        end if
        call jeveuo(lisnoe, 'L', jlinoe)
    end if
!
    if (lSelectDof) then
        chapil = sdpilo(1:14)//'.PLCR'
        call vtcreb(chapil, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
        call jeveuo(chapil(1:19)//'.VALE', 'E', jvale)
        call jeveuo(liscmp, 'L', jlicmp)
        call jelira(liscmp, 'LONMAX', ival=nddl)
!
        do iddl = 1, nddl
            nomcmp = zk8(jlicmp-1+iddl)
            do ino = 1, nbno
                if (lSelectDof) then
                    nume_node = zi(jlinoe-1+ino)
                    call nueqch('F', chapil, nume_node, nomcmp, numequ)
                end if

                if (lSelectDof) then
                    zr(jvale-1+numequ) = coef
                end if
            end do
        end do
    end if
!
    call jedetr(lisnoe)
    call jedetr(liscmp)

! --- CREATION SD REPERAGE DES DX/DY/DZ
    if (typpil .eq. 'LONG_ARC') then
        nb_dof_acti = 0
        selpil = sdpilo(1:14)//'.PLSL'
        call vtcreb(selpil, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
        call jeveuo(selpil(1:19)//'.VALE', 'E', vr=plsl)
        nddl = 3
        call wkvect(liscmp, 'V V K8', nddl, jlicmp)
        zk8(jlicmp-1+1) = 'DX'
        zk8(jlicmp-1+2) = 'DY'
        zk8(jlicmp-1+3) = 'DZ'
        do iddl = 1, nddl
            nomcmp = zk8(jlicmp-1+iddl)
            do ino = 1, nb_node_mesh
                nume_node = ino
                call nueqch(' ', selpil, nume_node, nomcmp, numequ)
                if (numequ .ne. 0) then
                    plsl(numequ) = 1.d0
                    nb_dof_acti = nb_dof_acti+1
                end if
            end do
        end do
        if (nb_dof_acti .eq. 0) then
            call utmess('F', 'MECANONLINE5_52')
        end if
    end if

! --- GESTION RECHERCHE LINEAIRE
    if (ds_algopara%l_line_search) then
        relmet = ds_algopara%line_search%method
        if (typpil .ne. 'DDL_IMPO') then
            if (relmet .ne. 'PILOTAGE') then
                call utmess('F', 'PILOTAGE_4')
            end if
        end if
    end if
!
    call jedetr(lisnoe)
    call jedetr(liscmp)
    call jedema()
end subroutine
