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

subroutine op0075()
    implicit none
!
!     OPERATEUR REST_GENE_PHYS
!
! ----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/elim75.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/harm75.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/regene.h"
#include "asterfort/rsadpa.h"
#include "asterfort/tran75.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!     ------------------------------------------------------------------
    character(len=8) :: k8bid, nomres, resin, mode, blanc8, param(3), val_param(3), ltch
    character(len=16) :: concep, nomcmd, typres, typrep, champ(8), typmat
    character(len=19) :: profno
    character(len=24) :: matgen, massgen, numgen, basemo
    aster_logical :: prsimp, l_err
    integer(kind=8) :: nbord, i, iord, lpain(3), lpaout(3), ibid, ir1
    integer(kind=8) :: j, j3refe, jrefn, naccab, nbcham
    integer(kind=8), pointer :: ordr(:) => null()
    character(len=24), pointer :: refa(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
    call infmaj()
!     ------------------------------------------------------------------
    k8bid = '        '
    blanc8 = '        '
    param(1) = 'MODELE'
    param(2) = 'CHAMPMAT'
    param(3) = 'CARAELEM'
!     -----------------------------------------------------------------
    call getres(nomres, typres, nomcmd)
!
!     --- PHASE DE TEST SUR LES CHAMPS A RESTITUER

!   on verifie que si ACCE_MONO_APPUI est présent TOUT_CHAM ou ACCE_ABSOLU est présent
    l_err = ASTER_FALSE
    call getvid(' ', 'ACCE_MONO_APPUI', scal=k8bid, nbret=naccab)
    if (naccab .gt. 0) l_err = ASTER_TRUE

    call getvtx(' ', 'NOM_CHAM', nbval=8, vect=champ, nbret=nbcham)
    if (nbcham .lt. 0) then
        call utmess('E', 'ALGORITH9_44')
    elseif (nbcham .eq. 0) then
!       TOUT_CHAM est présent
        l_err = ASTER_FALSE
    else
        do i = 1, nbcham
            do j = i+1, nbcham
                if (champ(i) .eq. champ(j)) then
                    call utmess('E', 'ALGORITH9_30')
                end if
            end do
            if (champ(i) .eq. 'ACCE_ABSOLU') then
                l_err = ASTER_FALSE
                if (naccab .eq. 0) then
                    call utmess('E', 'ALGORITH9_45')
                end if
            end if
        end do
    end if

    if (l_err) call utmess('F', 'ALGORITH9_38')
!
!     --- CREATION DU .REFN DU PROFIL :
    profno = '&&OP0075'//'.PROFC.NUME'
!
    call wkvect(profno(1:19)//'.REFN', 'V V K24', 5, jrefn)
    zk24(jrefn+1) = 'DEPL_R'
!
    call getvid(' ', 'RESU_GENE', scal=resin, nbret=ir1)
    call gettco(resin, concep)
!
!     --- INDICATEUR : 1) CALCUL CLASSIQUE AVEC UNE SIMPLE PROJECTION
!       -           OU 2) SANS MATRICE GENERALISEE (PROJ_MESU_MODAL)
    prsimp = .true.
!
    call dismoi('BASE_MODALE', resin, 'RESU_DYNA', repk=basemo, arret='C')
    call dismoi('NUME_DDL', resin, 'RESU_DYNA', repk=numgen, arret='C')
    call dismoi('REF_RIGI_PREM', resin, 'RESU_DYNA', repk=matgen, arret='C')
!
!   --- LE RESU_GENE NE VIENT PAS DE PROJ_MESU_MODAL
    if ((matgen(1:8) .ne. blanc8) .or. (numgen(1:8) .ne. blanc8)) then
        typrep = 'MODE_MECA'
        if (basemo(1:8) .ne. blanc8) call gettco(basemo, typrep)
!       --- LA BASE REFERENCEE DANS LE .REFD N'EST PAS UN MODE_MECA
        if (typrep(1:9) .ne. 'MODE_MECA') then
!           --- TRAITEMENT DE ELIM_LAGR
            if (typrep(1:16) .eq. 'MAILLAGE_SDASTER') then
                call gettco(matgen, typmat)
                ASSERT(typmat .eq. 'MATR_ASSE_ELIM_R')
                call dismoi('REF_MASS_PREM', resin, 'RESU_DYNA', repk=massgen)
                call gettco(massgen, typmat)
                ASSERT(typmat .eq. 'MATR_ASSE_ELIM_R')
            else
                prsimp = .false.
!             --- CHERCHER ALORS LE NUME_DDL_GENE POUR Y TROUVER DES INFOS
!               - SUR UN POTENTIEL MODE_GENE (CAS D'UNE DOUBLE-RESTITUTION)
                if (numgen(1:8) .eq. blanc8) then
!               --- PAS D'ENTREE DANS LE .REFD => CHERCHER DANS LA MATRICE K
                    call jeveuo(matgen(1:8)//'           .REFA', 'L', vk24=refa)
                    numgen = refa(2) (1:14)
                end if
                call jeveuo(numgen(1:14)//'.NUME.REFN', 'L', j3refe)
                call gettco(zk24(j3refe), typrep)
            end if
        end if
    end if
!
!
!
!     --- DYNAMIQUE TRANSITOIRE ---
    if (concep(1:9) .eq. 'TRAN_GENE') then
        if (prsimp) then
!         --- SIMPLE RESTITUTION => APPELER TRAN75 AVEC MODE=BLANC8
            call tran75(nomres, typres, resin, blanc8)
!
        else if (typrep(1:9) .eq. 'MODE_GENE') then
!         --- RECUPERER LA BASE MODALE POUR DOUBLE RESTITUTION
            call getvid(' ', 'MODE_MECA', scal=mode, nbret=ibid)
            if (ibid .eq. 0) then
                call utmess('F', 'ALGORITH9_48')
            end if
!
            call tran75(nomres, typres, resin, mode)
!
        else
!     REMARQUE 1: BLINDAGE
!         --- IMPOSSIBILE DE DETERMINER LE TYPE DE RESTITUTION A PARTIR
!           - DE LA SD_DYNA_GENE, LA SD A PROBABLEMENT ETE MAL DEFINIE
!           - A LA BASE. ON ARRETE LE CALCUL.
            ASSERT(.false.)
        end if
!
!     --- CALCUL MODAL AVEC ELIM_LAGR='OUI'
    else if (concep(1:9) .eq. 'MODE_MECA') then
        call elim75(nomres, resin, matgen, massgen)
!
!     --- CALCUL MODAL SANS SOUS-STRUCTURATION
    else if (concep(1:9) .eq. 'MODE_GENE') then
        call regene(nomres, resin, profno)
!
!     --- CALCUL HARMONIQUE
    else if (concep(1:9) .eq. 'HARM_GENE') then
        if (prsimp) then
!         --- SANS SOUS STRUCTURATION
            call harm75(nomres, typres, resin, blanc8)
!
        else if (typrep(1:9) .eq. 'MODE_GENE') then
!         --- AVEC SOUS STRUCTURATION, RECUPERER LA BASE MODALE
            call getvid(' ', 'MODE_MECA', scal=mode, nbret=ibid)
            if (ibid .eq. 0) then
                call utmess('F', 'ALGORITH9_48')
            end if
!
            call harm75(nomres, typres, resin, mode)
!
        else
!         --- BLINDAGE : VOIR REMARQUE 1
            ASSERT(.false.)
        end if
    end if
!
!     --- STOCKAGE DES RESULTATS
    call gettco(resin, concep)
    if (concep(1:9) .eq. 'MODE_GENE' .or. concep(1:9) .eq. 'MODE_MECA') then
        call jeveuo(nomres//'           .ORDR', 'L', vi=ordr)
        call jelira(nomres//'           .ORDR', 'LONUTI', nbord)
        do iord = 1, nbord
            call rsadpa(resin, 'L', 3, param, ordr(iord), &
                        0, tjv=lpain, styp=k8bid)
            call rsadpa(nomres, 'E', 3, param, ordr(iord), &
                        0, tjv=lpaout, styp=k8bid)
            do i = 1, 3
                zk8(lpaout(i)) = zk8(lpain(i))
            end do
        end do
    else
        call jeveuo(nomres//'           .ORDR', 'L', vi=ordr)
        call jelira(nomres//'           .ORDR', 'LONUTI', nbord)

        call dismoi('MODELE', basemo, 'RESULTAT', repk=val_param(1), arret='C')
        call dismoi('CHAM_MATER', basemo, 'RESULTAT', repk=val_param(2), arret='C')
        call dismoi('CARA_ELEM', basemo, 'RESULTAT', repk=val_param(3), arret='C')

        do i = 1, 3
            if (val_param(i) (1:6) .eq. '#AUCUN') val_param(i) = ' '
        end do

        do iord = 1, nbord
            call rsadpa(nomres, 'E', 3, param, ordr(iord), &
                        0, tjv=lpaout, styp=k8bid)
            do i = 1, 3
                zk8(lpaout(i)) = val_param(i)
            end do
        end do
    end if
!
    call jedema()
end subroutine
