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
subroutine tpsivp(p, sigmav)
    implicit none
#include "asterfort/assert.h"
#include "blas/dgemm.h"
#include "blas/dsymm.h"
    real(kind=8), dimension(:, :), intent(in) :: p
    real(kind=8), dimension(:), intent(inout) :: sigmav
!     ------------------------------------------------------------------
!     Appliquer un Changement de Base à un tenseur Sigma (contraintes ou
!                                                         déformation)
!     SIGMA <- P^T * SIGMA * P
!     ------------------------------------------------------------------
!     IN     P(3,3)     R   Matrice de passage de la base 1 à la base 2
!     INOUT  SIGMAV(6)  R   Tenseur de contraintes ou de déformations
!                           6 composantes : XX  YY  ZZ  XY  XZ  YZ
!                           en entrée: sigmav est exprimé dans la base 1
!                           en sortie: sigmav est exprimé dans la base 2
    !
!     Routine in place
!     ------------------------------------------------------------------
    integer(kind=8) :: ld
    real(kind=8) :: alpha, beta
    real(kind=8), dimension(3, 3) :: sigma, temp
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
!     ------------------------------------------------------------------
    ASSERT(size(p, 1) == 3)
!      ASSERT(size(sigmav)==6)
!     On décompacte le tenseur pour le stocker comme une matrice pleine symétrique 'L'
    sigma(:, :) = 0.d0
    sigma(1, 1) = sigmav(1)
    sigma(2, 2) = sigmav(2)
    sigma(3, 3) = sigmav(3)
    sigma(2, 1) = sigmav(4)
    sigma(3, 1) = sigmav(5)
    sigma(3, 2) = sigmav(6)
    !
! temp = sigma*P
    ld = 3
    alpha = 1.d0
    beta = 0.d0
    b_ldc = to_blas_int(ld)
    b_ldb = to_blas_int(ld)
    b_lda = to_blas_int(ld)
    b_m = to_blas_int(ld)
    b_n = to_blas_int(ld)
    call dsymm('L', 'L', b_m, b_n, alpha, &
               sigma, b_lda, p, b_ldb, beta, &
               temp, b_ldc)
! sigma <- P^T*temp
    b_ldc = to_blas_int(3)
    b_ldb = to_blas_int(3)
    b_lda = to_blas_int(3)
    b_m = to_blas_int(3)
    b_n = to_blas_int(3)
    b_k = to_blas_int(3)
    call dgemm('T', 'N', b_m, b_n, b_k, &
               alpha, p, b_lda, temp, b_ldb, &
               beta, sigma(1, 1), b_ldc)
! On re-compacte le tenseur
!xx
    sigmav(1) = sigma(1, 1)
! yy
    sigmav(2) = sigma(2, 2)
! zz
    sigmav(3) = sigma(3, 3)
! xy
    sigmav(4) = sigma(2, 1)
! xz
    sigmav(5) = sigma(3, 1)
! yz
    sigmav(6) = sigma(3, 2)
    !
end subroutine tpsivp
