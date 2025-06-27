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
subroutine pmathm(dimmat, dimdef, dimcon, dimuel, dsde, &
                  drds, ck, b, poids, work1, &
                  work2, matri)
    implicit none
#include "blas/dgemm.h"
    integer(kind=8) :: dimdef, dimcon, dimuel, dimmat
    real(kind=8) :: dsde(dimcon, dimdef), drds(dimdef, dimcon), poids
    real(kind=8) :: ck(dimdef), b(dimdef, dimuel)
    real(kind=8) :: work1(dimcon, dimuel), work2(dimdef, dimuel)
    real(kind=8) :: matri(dimmat, dimmat)
! ======================================================================
! --- BUT : PRODUIT DES MATRICES BT,C,DRDS,D,DSDE,F,B*POIDS ------------
! ---       CONTRIBUTION DU POINT D'INTEGRATION A DF -------------------
! ---       C,F,D SONT DIAGONALES --------------------------------------
! ======================================================================
    integer(kind=8) :: i, j
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
! ======================================================================
! --- ON FAIT LE CALCUL EN QUATRE FOIS ---------------------------------
! ======================================================================
!   WORK1 = DSDE x B
    b_ldc = to_blas_int(dimcon)
    b_ldb = to_blas_int(dimdef)
    b_lda = to_blas_int(dimcon)
    b_m = to_blas_int(dimcon)
    b_n = to_blas_int(dimuel)
    b_k = to_blas_int(dimdef)
    call dgemm('N', 'N', b_m, b_n, b_k, &
               1.d0, dsde, b_lda, b, b_ldb, &
               0.d0, work1, b_ldc)
!   WORK2 = DRDS x WORK1
    b_ldc = to_blas_int(dimdef)
    b_ldb = to_blas_int(dimcon)
    b_lda = to_blas_int(dimdef)
    b_m = to_blas_int(dimdef)
    b_n = to_blas_int(dimuel)
    b_k = to_blas_int(dimcon)
    call dgemm('N', 'N', b_m, b_n, b_k, &
               1.d0, drds, b_lda, work1, b_ldb, &
               0.d0, work2, b_ldc)
!   WORK2 = CK x WORK2
    do j = 1, dimuel
        do i = 1, dimdef
            work2(i, j) = ck(i)*work2(i, j)
        end do
    end do
!   MATRI = MATRI + POIDS x Bt x WORK2
    b_ldc = to_blas_int(dimmat)
    b_ldb = to_blas_int(dimdef)
    b_lda = to_blas_int(dimdef)
    b_m = to_blas_int(dimuel)
    b_n = to_blas_int(dimuel)
    b_k = to_blas_int(dimdef)
    call dgemm('T', 'N', b_m, b_n, b_k, &
               poids, b, b_lda, work2, b_ldb, &
               1.d0, matri, b_ldc)
! ======================================================================
end subroutine
