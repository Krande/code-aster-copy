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
subroutine op0079()
!
!
!  CALCUL PROJECTION SD_RESULTAT SUR BASE DE RITZ
!
!----------------------------------------------------------------------
!
    implicit none
!
!
!
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/copmod.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mdall2.h"
#include "asterfort/rrlds.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/trlds.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
!
    integer(kind=8) :: jsmde, nbmode, nbo, ii, iret, nbsym, idbase
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ibid, icod, iadref
    integer(kind=8) :: iddeeq, imod, ind, iord, isym
    integer(kind=8) :: jmod, n0, n1, n2, n4
    integer(kind=8) :: nbid, neq, tmod(1)
    real(kind=8) :: bid, ebid, pij
!-----------------------------------------------------------------------
    parameter(nbsym=3)
    character(len=1) :: typvec
    character(len=8) :: nomres, basemo, nomtyp, k8bid, res, numgen
    character(len=9) :: nosyin(nbsym)
    character(len=4) :: nosyou(nbsym), nosy
    character(len=14) :: nu, numdd1, numdd2
    character(len=16) :: typres, nomcom, typbas, matri2
    character(len=19) :: nochno
    character(len=24) :: matric, deeq
    complex(kind=8) :: cbid
    real(kind=8), pointer :: matrnorm(:) => null()
    real(kind=8), pointer :: vectass1(:) => null()
    real(kind=8), pointer :: vectass2(:) => null()
    real(kind=8), pointer :: vectasse(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: ordr(:) => null()
    integer(kind=8), pointer :: nequ(:) => null()
    blas_int :: b_incx, b_incy, b_n
    data nosyin/'DEPL', 'VITE', 'ACCE'/
    data nosyou/'DEPL', 'VITE', 'ACCE'/
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
! REMARQUE : ACTUELLEMENT, SEULS LES CHAMPS DEPL, VITE ET ACCE SONT
! TRAITES. IL CONVIENDRAIT NORMALEMENT DE TRAITER LES FORC_NODAL, MAIS
! POUR CELA, IL FAUT CREER UN CHAMP EQUIVALENT DANS LA SD_TRAN_GENE.
! OBJECTIF : POUVOIR DEFINIR UN CHARGEMENT GENERALISE POUR DTM AVEC UN
! TYPAGE CORRECT POUR LE CHARGEMENT. ACTUELLEMENT, LE CHARGEMENT
! APPLIQUE AVEC EXCIT_RESU EST UNE SD RESULTAT AVEC DES CHAMPS DEPL,
! DANS DYNA_VIBRA.
!
    call jemarq()
    call infmaj()
!
! --- RECUPERATION DES ARGUMENTS DE LA COMMANDE
!
    call getres(nomres, typres, nomcom)
    call getvid(' ', 'NUME_DDL_GENE', scal=numgen, nbret=n0)
    call getvid(' ', 'RESU', scal=res, nbret=n1)
! LE CAS RESU_GENE N'EST PAS ACTIVE POUR LE MOMENT
!      CALL GETVID(' ','RESU_GENE',0,IARG,1,RES,N3)
    call getvid(' ', 'BASE', scal=basemo, nbret=n4)
    call getvtx(' ', 'TYPE_VECT', scal=nomtyp, nbret=n2)
    call gettco(basemo, typbas)
!
! --- RECUPERATION DU NB DE MODES
!
    call rsorac(basemo, 'LONUTI', ibid, bid, k8bid, &
                cbid, ebid, 'ABSOLU', tmod, 1, &
                nbid)
    nbmode = tmod(1)
!
!
    call jeveuo(numgen//'      .SMOS.SMDE', 'L', jsmde)
!
!
! --- RECUPERATION DU NOMBRE DE NUME_ORDRE DE LA SD_RESU
!
    call rsorac(res, 'LONUTI', ibid, bid, k8bid, &
                cbid, ebid, 'ABSOLU', tmod, 1, &
                nbid)
    nbo = tmod(1)
    call jeveuo(res//'           .ORDR', 'L', vi=ordr)
!
!
! --- VERIFICATION DE LA CONFORMITE DES NUMEROTATIONS
!     DES MODES ET DU VECTEUR ASSEMBLE
!     ON RECUPERE LES NUME_DDL DANS LES REFD DES DEUX SD
!     SI ELLES SONT ABSENTES, ON ESSAYE AVEC LES MATRICES
!
!
!
!
    if (typbas(1:9) .eq. 'MODE_MECA') then
        call dismoi('NUME_DDL', res, 'RESU_DYNA', repk=nu)
        if (nu(1:1) .ne. ' ') then
            numdd1 = nu
        else
            call dismoi('REF_RIGI_PREM', res, 'RESU_DYNA', repk=matric, arret='C')
            call exisd('MATR_ASSE', matric, iret)
            if (iret .ne. 0) then
                call dismoi('NOM_NUME_DDL', matric, 'MATR_ASSE', repk=nu)
                numdd1 = nu
            end if
            if (iret .eq. 0) then
                call utmess('F', 'ALGORITH17_8', sk=res)
            end if
        end if
        call dismoi('NUME_DDL', basemo, 'RESU_DYNA', repk=nu)
        if (nu(1:1) .ne. ' ') then
            numdd2 = nu
        else
            call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric, arret='C')
            call exisd('MATR_ASSE', matric, iret)
            if (iret .ne. 0) then
                call dismoi('NOM_NUME_DDL', matric, 'MATR_ASSE', repk=nu)
                numdd2 = nu
            end if
            if (iret .eq. 0) then
                call utmess('F', 'ALGORITH17_8', sk=basemo)
            end if
        end if
!
    else if (typbas(1:9) .eq. 'MODE_GENE') then
        call dismoi('NUME_DDL', res, 'RESU_DYNA', repk=nu)
        call dismoi('REF_RIGI_PREM', basemo, 'RESU_DYNA', repk=matric)
        matri2 = matric(1:16)
        call jeveuo(matri2//'   .REFA', 'L', vk24=refa)
        numdd2 = refa(2) (1:14)
    end if
!
    if (numdd1 .ne. numdd2) then
        call utmess('I', 'ALGORITH9_41')
    end if
!
! --- RECUPERATION DU NOMBRE D'EQUATIONS DU SYSTEME PHYSIQUE
!
    if ((typbas(1:9) .eq. 'MODE_MECA')) then
        call dismoi('NB_EQUA', numdd1, 'NUME_DDL', repi=neq)
    else if (typbas(1:9) .eq. 'MODE_GENE') then
        call jeveuo(numdd1//'.NUME.NEQU', 'L', vi=nequ)
        neq = nequ(1)
    end if
!
    deeq = nu//'.NUME.DEEQ'
    call jeveuo(deeq, 'L', iddeeq)
!
! --- INITIALISATION DE LA SD_RESULTAT
!
    call mdall2(nomres, basemo, res, nbo, nbmode)
!
! --- RECUPERE LA BASE MODALE SOUS LA FORME D'UN VECT NBMODE*NEQ
!
    call wkvect('&&OP0072.BASEMO', 'V V R', nbmode*neq, idbase)
    call copmod(basemo, numer=nu, bmodr=zr(idbase))
!
! --- BOUCLE SUR LES NUM_ORDR ET LES NOMSY DE LA SD_RESULTAT
!     ATTENTION : ON NE TRAITE QUE LES NOMSY STOCKABLES DANS
!     UN TRAN_GENE : DEPL, ACCE, VITE
!
    do isym = 1, nbsym
!
        do iord = 1, nbo
!
            nosy = nosyou(isym)
!
! --- RECUP DU CHAMP DE LA SDIN CORRESPONDANT AU NUME_ORDR ET ISYM
            call rsexch(' ', res, nosyin(isym), ordr(iord), nochno, &
                        iret)
            if (iret .ne. 0) goto 40
            call jeveuo(nochno//'.VALE', 'L', vr=vale)
            call jeveuo(nochno//'.REFE', 'L', iadref)
            call jelira(nochno//'.VALE', 'TYPE', cval=typvec)
! --- LE CAS COMPLEXE (SD HARMONIQUES) N'EST PAS TRAITE
            if (typvec .eq. 'C') then
                call utmess('F', 'ALGORITH17_19')
            end if
!
! --- INDICE DE STOCKAGE
            call jeveuo(nomres//'           .'//nosy, 'E', ii)
!
            if (nomtyp(1:4) .eq. 'FORC') then
!
! --- PROJECTION D UN VECTEUR DE TYPE FORCE
!
                AS_ALLOCATE(vr=vectasse, size=neq)
                do imod = 1, nbmode
!
!
! --------- RECOPIE DU IEME MODE DANS UN VECTEUR TEMP
!
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(idbase+(imod-1)*neq), b_incx, vectasse, b_incy)
!
! ------- MISE A ZERO DES DDLS DE LAGRANGE
!
                    call zerlag(neq, zi(iddeeq), vectr=vectasse)
!
! ------- PRODUIT SCALAIRE VECTASS * MODE
!
                    ind = ii-1+(iord-1)*nbmode+imod
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    zr(ind) = ddot(b_n, vectasse, b_incx, vale, b_incy)
!
! ------- LIBERATION DU VECTEUR TEMP
                end do
                AS_DEALLOCATE(vr=vectasse)
            else
!
! --- PROJECTION D UN VECTEUR DE TYPE DEPL OU VITE
!
                AS_ALLOCATE(vr=vectass1, size=neq)
                AS_ALLOCATE(vr=vectass2, size=neq)
                AS_ALLOCATE(vr=matrnorm, size=nbmode*nbmode)
!
! ----- CALCUL DE TMODE*MODE
!
                do imod = 1, nbmode
!
! ----- RECOPIE DU IEME MODE
!
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(idbase+(imod-1)*neq), b_incx, vectass1, b_incy)
!
! ------- MISE A ZERO DES DDLS DE LAGRANGE
!
                    call zerlag(neq, zi(iddeeq), vectr=vectass1)
!
!-------- PRODUIT SCALAIRE MODE(IMOD)*MODE(JMOD)
!
                    do jmod = imod, nbmode
!
! ------- RECOPIE DU JEME MODE
!
                        b_n = to_blas_int(neq)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        call dcopy(b_n, zr(idbase+(jmod-1)*neq), b_incx, vectass2, b_incy)
! --------- MISE A ZERO DES DDLS DE LAGRANGE
!
                        call zerlag(neq, zi(iddeeq), vectr=vectass2)
!
! --------- PRODUIT SCALAIRE MODE(IMOD)*MODE(JMOD)
!
                        b_n = to_blas_int(neq)
                        b_incx = to_blas_int(1)
                        b_incy = to_blas_int(1)
                        pij = ddot(b_n, vectass1, b_incx, vectass2, b_incy)
                        matrnorm(1+imod+(jmod-1)*nbmode-1) = pij
                        matrnorm(1+jmod+(imod-1)*nbmode-1) = pij
                    end do
                end do
!
! ----- CALCUL DE LA PROJECTION
!
                do imod = 1, nbmode
!
! ------- RECOPIE DU IEME MODE
!
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(idbase+(imod-1)*neq), b_incx, vectass1, b_incy)
!
! ------- MISE A ZERO DES DDLS DE LAGRANGE
!
                    call zerlag(neq, zi(iddeeq), vectr=vectass1)
!
! ------- PRODUIT SCALAIRE VECTASS * MODE
!
                    b_n = to_blas_int(neq)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    vectass2(imod) = ddot(b_n, vectass1, b_incx, vale, b_incy)
!
                end do
!
!
! ----- FACTORISATION ET RESOLUTION SYSTEME
!
                ind = ii+(iord-1)*nbmode
                call trlds(matrnorm, nbmode, nbmode, icod)
                if (icod .ne. 0) then
                    call utmess('F', 'ALGORITH9_42')
                end if
                call rrlds(matrnorm, nbmode, nbmode, vectass2, 1)
                b_n = to_blas_int(nbmode)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, vectass2, b_incx, zr(ind), b_incy)
                AS_DEALLOCATE(vr=vectass1)
                AS_DEALLOCATE(vr=vectass2)
                AS_DEALLOCATE(vr=matrnorm)
                if (typvec .eq. 'C') call jedetr('&&OP0079.VECTASC2')
            end if
!
        end do
40      continue
    end do
!
    call jedema()
end subroutine
