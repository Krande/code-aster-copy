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
subroutine vecind(mat, lvec, nbl, nbc, force, &
                  nindep)
    implicit none
!
!--------------------------------------------------------------------C
!  M. CORUS     DATE 09/06/11
!-----------------------------------------------------------------------
!  BUT : CONSTRUCTION D'UNE BASE A PARTIR DE VECTEURS QUELCONQUES :
!         - SELECTION D'UNE FAMILLE LIBRE
!         - ORTHONORMALISATION DE LA FAMILLE LIBRE
!
!
!  MAT     /I/   : NOM K19 DE LA MATRICE POUR CONSTRUIRE LA NORME
!  LVEC    /I-O/ : POINTEUR DE LA FAMILLE DE VECTEURS
!  NBL     /I/   : NOMBRE DE LIGNE DE CHAQUE VECTEUR
!  NBC     /I-O/ : NOMBRE DE VECTEURS
!  FORCE   /I/   : FORCE LA NORMALISATION SI VAUT 1
!  NINDEP  /O/   : NOMBRE DE VECTEUR INDEPENDANT EN SORTIE
!
!  NOTE : LA TAILLE DE LA MATRICE ASSOCIE A LVEC EST INCHANGEE EN SORTIE
!         MAIS LES DERNIERES COLONNES SONT MISES A ZERO
!-----------------------------------------------------------------------
!
!
!
!
!
#include "jeveux.h"
#include "asterc/matfpe.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lceqvn.h"
#include "asterfort/mrmult.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/daxpy.h"
#include "blas/ddot.h"
#include "blas/dgemm.h"
#include "blas/dgesvd.h"
    integer(kind=8) :: lvec, nbl, nbc, nindep, lwork, lmat, ltrav1
    integer(kind=8) :: i1, k1, l1, iret, lcopy, force, indnz
    integer(kind=4) :: info
    real(kind=8) :: swork(1), norme, sqrt, rij
    character(len=8) :: ortho
    character(len=19) :: mat, nume
    real(kind=8), pointer :: mat_svd_work(:) => null()
    real(kind=8), pointer :: new_stat(:) => null()
    real(kind=8), pointer :: trav2_u(:) => null()
    real(kind=8), pointer :: trav3_v(:) => null()
    integer(kind=8), pointer :: vec_ind_nz(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    blas_int :: b_incx, b_incy, b_n
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m
    blas_int :: b_ldu, b_ldvt, b_lwork
!
    ortho = ' '
    iret = 0
    AS_ALLOCATE(vr=new_stat, size=nbc*nbc)
    call wkvect('&&VECIND.TRAV1', 'V V R', nbl, ltrav1)
    AS_ALLOCATE(vi=vec_ind_nz, size=nbc)
    indnz = 0
    if (mat .ne. ' ') then
        call jeveuo(mat//'.&INT', 'L', lmat)
        call dismoi('NOM_NUME_DDL', mat(1:8), 'MATR_ASSE', repk=nume)
        call jeveuo(nume(1:8)//'      .NUME.DEEQ', 'L', vi=deeq)
    end if
!
!-- NORMER LES MODES DANS L2 AVANT DE CONSTRUIRE LA MATRICE
!-- POUR PLUS DE ROBUSTESSE
!
    call wkvect('&&VECIND.VECTEURS_COPIES', 'V V R', nbl*nbc, lcopy)
    call wkvect('&&VECIND.VECTEURS_TEMP', 'V V R', nbl, ltrav1)
!
    do i1 = 1, nbc
        if (mat .ne. ' ') then
            call zerlag(nbl, deeq, vectr=zr(lvec+nbl*(i1-1)))
            call mrmult('ZERO', lmat, zr(lvec+nbl*(i1-1)), zr(ltrav1), 1, &
                        .true._1)
            call zerlag(nbl, deeq, vectr=zr(ltrav1))
            b_n = to_blas_int(nbl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            norme = ddot(b_n, zr(ltrav1), b_incx, zr(lvec+nbl*(i1-1)), b_incy)
!
        else
            b_n = to_blas_int(nbl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            norme = ddot(b_n, zr(lvec+nbl*(i1-1)), b_incx, zr(lvec+nbl*(i1-1)), b_incy)
        end if
        norme = sqrt(norme)
        if (norme .gt. 1.d-16) then
            b_n = to_blas_int(nbl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1/norme, zr(lvec+nbl*(i1-1)), b_incx, zr(lcopy+nbl*(i1-1)), &
                       b_incy)
!        ELSE
!          CALL DAXPY(NBL,0.D0,ZR(LVEC+NBL*(I1-1)),1,
!     &               ZR(LCOPY+NBL*(I1-1)),1)
        end if
    end do
!
    do l1 = 1, nbc
        if (mat .ne. ' ') then
            call mrmult('ZERO', lmat, zr(lcopy+nbl*(l1-1)), zr(ltrav1), 1, &
                        .true._1)
        else
            call lceqvn(nbl, zr(lcopy+nbl*(l1-1)), zr(ltrav1))
        end if
        b_n = to_blas_int(nbl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        norme = ddot(b_n, zr(ltrav1), b_incx, zr(lcopy+nbl*(l1-1)), b_incy)
        new_stat(1+(l1-1)*(nbc+1)) = norme
        do k1 = l1+1, nbc
            b_n = to_blas_int(nbl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            rij = ddot(b_n, zr(ltrav1), b_incx, zr(lcopy+nbl*(k1-1)), b_incy)
            new_stat(1+(l1-1)*nbc+k1-1) = rij
            new_stat(1+(k1-1)*nbc+l1-1) = rij
        end do
    end do
!
    if (force .ne. 1) then
!-- UTILISE APRES LE GRAM SCHMIDT DANS ORTH99, POUR
!-- ELIMINER LES VECTEURS NON INDEPENDANTS
!
        do l1 = 1, nbc
            norme = new_stat(1+(l1-1)*(nbc+1))
!
            if (norme .gt. 1.d-16) then
                do k1 = l1+1, nbc
                    rij = abs(new_stat(1+(l1-1)*nbc+k1-1))
                    rij = rij/norme
                    if (rij .gt. 1.d-8) then
                        write (6, *) ' ... ANNULATION DU VECTEUR ', k1
                        do i1 = 1, nbl
                            zr(lvec+((k1-1)*nbl)+i1-1) = 0.d0
                        end do
                        do i1 = 1, nbc
                            new_stat(1+((k1-1)*nbc)+i1-1) = 0.d0
                            new_stat(1+((i1-1)*nbc)+k1-1) = 0.d0
                        end do
                    end if
                end do
            end if
!
        end do
!
        call getvtx('  ', 'ORTHO', iocc=1, nbval=8, vect=ortho, &
                    nbret=iret)
        if ((iret .eq. 1) .and. (ortho .eq. 'OUI')) then
!-- SELECTION DES VECTEURS NON NULS POUR REMPLIR LA BASE
            do i1 = 1, nbc
                if (new_stat(1+(i1-1)*(nbc+1)) .gt. 1d-10) then
                    vec_ind_nz(indnz+1) = i1
                    indnz = indnz+1
                end if
            end do
!
            do i1 = 1, indnz
                l1 = vec_ind_nz(i1)
                if (i1 .ne. l1) then
                    call lceqvn(nbl, zr(lvec+nbl*(l1-1)), zr(lvec+nbl*(i1-1)))
                end if
            end do
            nbc = indnz
        end if
    else
!
!-- ALLOCATION DES MATRICES DE TRAVAIL TEMPORAIRES
        call wkvect('&&VECIND.TRAV1_S', 'V V R', nbc, ltrav1)
        AS_ALLOCATE(vr=trav2_u, size=nbc*nbc)
        AS_ALLOCATE(vr=trav3_v, size=nbc*nbc)
!
!-- DESACTIVATION DU TEST FPE
        call matfpe(-1)
!
        b_ldvt = to_blas_int(nbc)
        b_ldu = to_blas_int(nbc)
        b_lda = to_blas_int(nbc)
        b_m = to_blas_int(nbc)
        b_n = to_blas_int(nbc)
        b_lwork = to_blas_int(-1)
        call dgesvd('A', 'N', b_m, b_n, new_stat, &
                    b_lda, zr(ltrav1), trav2_u, b_ldu, trav3_v, &
                    b_ldvt, swork, b_lwork, info)
        lwork = int(swork(1))
        AS_ALLOCATE(vr=mat_svd_work, size=lwork)
!
        b_ldvt = to_blas_int(nbc)
        b_ldu = to_blas_int(nbc)
        b_lda = to_blas_int(nbc)
        b_m = to_blas_int(nbc)
        b_n = to_blas_int(nbc)
        b_lwork = to_blas_int(lwork)
        call dgesvd('A', 'N', b_m, b_n, new_stat, &
                    b_lda, zr(ltrav1), trav2_u, b_ldu, trav3_v, &
                    b_ldvt, mat_svd_work, b_lwork, info)
!
        nindep = 0
        norme = (nbc+0.d0)*zr(ltrav1)*1.d-16
        do k1 = 1, nbc
            if (zr(ltrav1+k1-1) .gt. norme) nindep = nindep+1
        end do
!
        call wkvect('&&VECIND.MODE_INTF_DEPL', 'V V R', nbl*nbc, lmat)
!
        b_ldc = to_blas_int(nbl)
        b_ldb = to_blas_int(nbc)
        b_lda = to_blas_int(nbl)
        b_m = to_blas_int(nbl)
        b_n = to_blas_int(nindep)
        b_k = to_blas_int(nbc)
        call dgemm('N', 'N', b_m, b_n, b_k, &
                   1.d0, zr(lcopy), b_lda, trav2_u, b_ldb, &
                   0.d0, zr(lvec), b_ldc)
!
!-- INUTILE D'ANNULER DES VECTEURS QUI NE SERVIRONT NUL PART...
!        DO 540 I1=1,NBL*(NBC-NINDEP)
!          ZR(LVEC+(NINDEP*NBL)+I1-1)=0.D0
! 540    CONTINUE
!
        call matfpe(1)
!
        AS_DEALLOCATE(vr=mat_svd_work)
        call jedetr('&&VECIND.TRAV1_S')
        AS_DEALLOCATE(vr=trav2_u)
        AS_DEALLOCATE(vr=trav3_v)
        call jedetr('&&VECIND.MODE_INTF_DEPL')
!
    end if
!
!
!
    AS_DEALLOCATE(vr=new_stat)
    call jedetr('&&VECIND.TRAV1')
    call jedetr('&&VECIND.VECTEURS_COPIES')
    call jedetr('&&VECIND.VECTEURS_TEMP')
    AS_DEALLOCATE(vi=vec_ind_nz)
!
end subroutine
