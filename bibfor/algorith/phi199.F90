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

subroutine phi199(model, mate, mateco, ma, nu, num, &
                  nbmode, solvez, indice, tabad)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/calflu.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/majou.h"
#include "asterfort/prstoc.h"
#include "asterfort/pteddl.h"
#include "asterfort/reliem.h"
#include "asterfort/resoud.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/rsvpar.h"
#include "asterfort/tabcor.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    integer(kind=8) :: nbmode, indice, tabad(*)
    character(len=8) :: ma
    character(len=14) :: nu, num
    character(len=*) :: mate, mateco, solvez
!
!
! CALCULS DES CONDITIONS AUX LIMITES POUR LA DETERMINATION
! DES POTENTIELS FLUCTUANTS POUR LA FORCE AJOUTEE,
! IN : K* : MODEL : TYPE DE MODELISATION FLUIDE
! IN : K* : MATE : MATERIAU FLUIDE
! IN : K* : PHIBAR : NOM DU POTENTIEL PERMANENT
! IN : K* : MA : NOM DE LA MATRICE DE RAIDEUR FLUIDE
! IN : K* : NU : NUMEROTATION DES DDLS ASSOCIES AU FLUIDE
! IN : K* : NUM : NUMEROTATION DES DDLS ASSOCIES A L'INTERFACE
!           POTENTIELS FLUCTUANTS : 1 : MASSE AJOUTEE
!                                 : 2 : AMORTISSEMENT ET RAIDEUR
! IN : K* : SOLVEZ : METHODE DE RESOLUTION 'MULT_FRONT','LDLT' OU 'GCPC'
!---------------------------------------------------------------------
    integer(kind=8) :: ibid, nbvale, nbrefe, iret, nbno, id, ier
    integer(kind=8) :: ilires, jref, neq, nbd, nbdir, i, jvec, in, nbsta
    integer(kind=8) :: iphi1, n3, n1, icor(2), n2, ndble, iordr, nbtrou, tmod(1)
    real(kind=8) :: rbid, xnorm, xd, depl(6), epsi
    complex(kind=8) :: c16b, cbid
    character(len=2) :: model
    character(len=8) :: k8bid, modmec, mailla, maflui, tabcmp(6), crit
    character(len=8) :: moflui, moint, typmcl(2), modsta
    character(len=14) :: nume
    character(len=16) :: acces, motcle(2)
    character(len=19) :: vecso1, vesto1, maprec, solveu, chsol, chamno
    character(len=24) :: nomcha, nocham, criter
    character(len=24) :: valk(3)
    integer(kind=8), pointer :: ddl(:) => null()
    real(kind=8), pointer :: mst(:) => null()
    character(len=8), pointer :: noeud(:) => null()
!
    data maprec/'&&OP0199.MAPREC'/
    data chsol/'&&OP0199.SOLUTION'/
    data tabcmp/'DX', 'DY', 'DZ', 'DRX', 'DRY', 'DRZ'/
    data ndble/0/
! -----------------------------------------------------------------
!
    call jemarq()
    epsi = r8prem()
    ier = 0
    solveu = solvez
    criter = '&&RESGRA_GCPC'
    indice = 0
!
    call getvid(' ', 'MODELE_FLUIDE', scal=moflui, nbret=n1)
    call getvid(' ', 'MODELE_INTERFACE', scal=moint, nbret=n2)
    call getvid(' ', 'MODE_MECA', scal=modmec, nbret=n3)
!
! --- TEST POUR DETERMINER SI FLUIDE ET STRUCTURE S APPUIENT SUR
!     DES MAILLAGES COMMUNS
!
    if (n3 .gt. 0) then
        call rsorac(modmec, 'LONUTI', 0, rbid, k8bid, &
                    cbid, rbid, 'ABSOLU', tmod, 1, &
                    ibid)
        nbmode = tmod(1)
        call rsexch('F', modmec, 'DEPL', 1, nomcha, &
                    iret)
        call dismoi('NOM_MAILLA', nomcha, 'CHAM_NO', repk=mailla)
        call dismoi('NOM_MAILLA', moint, 'MODELE', repk=maflui)
        if (maflui .ne. mailla) then
            call tabcor(model, mate, mateco, mailla, maflui, moint, &
                        num, ndble, icor)
            call majou(model, modmec, solveu, num, nu, &
                       ma, mate, mateco, moint, ndble, icor, &
                       tabad)
            indice = 1
        end if
    end if
!
!=====================================================================
!---------------- ALTERNATIVE MODE_MECA OU---------
!-----------------------------MODELE-GENE--------------------
!=====================================================================
! DANS LE CAS OU ON N A PAS CALCUL DE MASSE AJOUTEE SUR UN
! MAILLAGE SQUELETTE
!
    if ((n3 .gt. 0) .and. (indice .ne. 1)) then
!
!
!----- -RECUPERATION DU NB DE MODES DU CONCEPT MODE_MECA
!
        call rsorac(modmec, 'LONUTI', 0, rbid, k8bid, &
                    cbid, rbid, 'ABSOLU', tmod, 1, &
                    ibid)
        nbmode = tmod(1)
!
        call wkvect('&&OP0199.PHI1', 'V V K24', 1, iphi1)
!
!======================================================================
! BOUCLE SUR LE NOMBRE DE MODES: CALCUL DU FLUX FLUIDE MODAL
!======================================================================
        ilires = 0
        nomcha = '&&PHI199.CHAMREF'
        call rsexch('F', modmec, 'DEPL', 1, nocham, &
                    iret)
        nocham = nocham(1:19)//'.REFE'
        call jeveuo(nocham, 'L', jref)
        nume = zk24(jref+1) (1:14)
        call vtcreb(nomcha, 'V', 'R', &
                    nume_ddlz=nume, &
                    nb_equa_outz=neq)
!
! --- QUELLE EST LA DIRECTION ?
!
        call getvr8(' ', 'DIRECTION', nbval=0, nbret=nbd)
        nbdir = -nbd
        call getvr8(' ', 'DIRECTION', nbval=nbdir, vect=depl, nbret=nbd)
!
!     --- ON NORMALISE LE VECTEUR ---
        xnorm = 0.d0
        do i = 1, nbdir
            xnorm = xnorm+depl(i)*depl(i)
        end do
        xnorm = sqrt(xnorm)
        if (xnorm .lt. 0.d0) then
            call utmess('F', 'ALGORITH9_81')
        end if
        do i = 1, nbdir
            depl(i) = depl(i)/xnorm
        end do
!
        call jeveuo(nomcha(1:19)//'.VALE', 'E', jvec)
        AS_ALLOCATE(vi=ddl, size=neq*nbdir)
        call pteddl('NUME_DDL', nume, nbdir, tabcmp, neq, &
                    tabl_equa=ddl)
!
        do in = 0, neq-1
            zr(jvec+in) = 0.d0
        end do
!
!     --- ON RECUPERE LES MODES STATIQUES ---
!
        call getvid(' ', 'MODE_STAT', scal=modsta, nbret=nbsta)
        if (nbsta .eq. 0) goto 41
!
!     --- ON RECUPERE LES POINTS D'ANCRAGE ---
!
        motcle(1) = 'NOEUD'
        motcle(2) = 'GROUP_NO'
        typmcl(1) = 'NOEUD'
        typmcl(2) = 'GROUP_NO'
        call reliem(' ', mailla, 'NO_NOEUD', ' ', 1, &
                    2, motcle, typmcl, '&&PHI199.NOEUD', nbno)
        call jeveuo('&&PHI199.NOEUD', 'L', vk8=noeud)
!
!     --- ON BOUCLE SUR LES NOEUDS ---
!
        do id = 1, nbdir
            xd = depl(id)
            if (abs(xd) .gt. epsi) then
                do in = 1, nbno
                    acces(1:8) = noeud(in)
                    acces(9:16) = tabcmp(id)
!
!              --- ON RECUPERE LE MODE STATIQUE ASSOCIE AU NOEUD ---
                    call rsorac(modsta, 'NOEUD_CMP', ibid, rbid, acces, &
                                c16b, epsi, crit, tmod, 1, &
                                nbtrou)
                    iordr = tmod(1)
                    if (nbtrou .ne. 1) then
                        ier = ier+1
                        valk(1) = acces(1:8)
                        valk(2) = acces(9:16)
                        call utmess('E', 'ALGELINE4_61', nk=2, valk=valk)
                        goto 26
                    end if
                    call rsvpar(modsta, iordr, 'TYPE_DEFO', ibid, rbid, &
                                'DEPL_IMPO', iret)
                    if (iret .ne. 100) then
                        ier = ier+1
                        valk(1) = 'DDL_IMPO'
                        valk(2) = acces(1:8)
                        valk(3) = acces(9:16)
                        call utmess('E', 'ALGELINE4_62', nk=3, valk=valk)
                        goto 26
                    end if
                    call rsexch('F', modsta, 'DEPL', iordr, chamno, &
                                iret)
                    call jeveuo(chamno//'.VALE', 'L', vr=mst)
!
                    do i = 0, neq-1
                        zr(jvec+i) = zr(jvec+i)-ddl(1+(id-1)*neq+i)*xd*mst(1+i)
                    end do
                    call jelibe(chamno//'.VALE')
26                  continue
                end do
            end if
        end do
        if (ier .ne. 0) then
            call utmess('F', 'ALGORITH5_24')
        end if
!
        goto 42
!
41      continue
        do i = 1, nbdir
            do in = 0, neq-1
                zr(jvec+in) = zr(jvec+in)-ddl(1+(i-1)*neq+in)*depl(i)
            end do
        end do
42      continue
!
        nomcha = nomcha(1:19)
        vecso1 = '&&OP0199.VECSOL1'
!
        call calflu(nomcha, moflui, mate, mateco, nu, vecso1, &
                    nbrefe, nbvale, 'R')
!
        ilires = ilires+1
!
!------------- RESOLUTION  DU LAPLACIEN EN 2D-----------------------
!
        call resoud(ma, maprec, solveu, ' ', 0, &
                    vecso1, chsol, 'V', [0.d0], [cbid], &
                    criter, .true._1, 0, iret)
        call jedupc('V', chsol(1:19), 1, 'V', vecso1(1:19), &
                    .false._1)
        call detrsd('CHAMP_GD', chsol)
!
!------------ CREATION DU VECTEUR PRESSION MODAL-------------------
!
!- FORMATION DU TABLEAU CONTENANT LA PRESSION POUR CHAQUE MODE-----
!
!------------------------------------------------------------------
        vesto1 = '&&OP0199.VEST1'
        call prstoc(vecso1, vesto1, ilires, ilires, iphi1, &
                    nbvale, nbrefe)
!
    end if
!
! --- MENAGE
    call detrsd('CHAM_NO', '&&PHI199.CHAMREF')
    AS_DEALLOCATE(vi=ddl)
    call jedetr('&&PHI199.NOEUD')
!
    call jeexin(criter(1:19)//'.CRTI', iret)
    if (iret .ne. 0) then
        call jedetr(criter(1:19)//'.CRTI')
        call jedetr(criter(1:19)//'.CRTR')
        call jedetr(criter(1:19)//'.CRDE')
    end if
!----------------------------------------------------------------
    call jedema()
end subroutine
