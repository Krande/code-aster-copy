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
!
module HHO_algebra_module
!
    use HHO_matrix_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
#include "blas/dsymv.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - linear algebra routine
!
! --------------------------------------------------------------------------------------------------
!

    public :: hho_dgemm_NN, hho_dgemm_TN, hho_dgemv_N, hho_dgemv_T
    public :: hho_dsymv_U
!
contains
!---------------------------------------------------------------------------------------------------
! -- member function for HHO_matrix type
!---------------------------------------------------------------------------------------------------
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_dgemm_NN(alpha, matA, matB, beta, matC)
!
        implicit none
!
        type(HHO_matrix), intent(in) :: matA, matB
        type(HHO_matrix), intent(inout) :: matC
        real(kind=8), intent(in) :: alpha, beta
!
! --------------------------------------------------------------------------------------------------
!
!   DGEMM blas routine for hho_matrix
! --------------------------------------------------------------------------------------------------
!
        blas_int :: b_k, b_n, b_m, b_lda, b_ldb, b_ldc
!
        b_lda = to_blas_int(matA%max_nrows)
        b_ldb = to_blas_int(matB%max_nrows)
        b_ldc = to_blas_int(matC%max_nrows)
        b_m = to_blas_int(matA%nrows)
        b_n = to_blas_int(matB%ncols)
        b_k = to_blas_int(matA%ncols)
!
        ASSERT(matB%nrows == matA%ncols)
        ASSERT(matA%nrows == matC%nrows)
        ASSERT(matB%ncols == matC%ncols)
!
        call dgemm('N', 'N', b_m, b_n, b_k, &
                   alpha, matA%m, b_lda, matB%m, b_ldb, &
                   beta, matC%m, b_ldc)
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_dgemm_TN(alpha, matA, matB, beta, matC)
!
        implicit none
!
        type(HHO_matrix), intent(in) :: matA, matB
        type(HHO_matrix), intent(inout) :: matC
        real(kind=8), intent(in) :: alpha, beta
!
! --------------------------------------------------------------------------------------------------
!
!   DGEMM blas routine for hho_matrix
! --------------------------------------------------------------------------------------------------
!
        blas_int :: b_k, b_n, b_m, b_lda, b_ldb, b_ldc
!
        b_lda = to_blas_int(matA%max_nrows)
        b_ldb = to_blas_int(matB%max_nrows)
        b_ldc = to_blas_int(matC%max_nrows)
        b_m = to_blas_int(matA%ncols)
        b_n = to_blas_int(matB%ncols)
        b_k = to_blas_int(matA%nrows)
!
        ASSERT(matB%nrows == matA%nrows)
        ASSERT(matA%ncols == matC%nrows)
        ASSERT(matB%ncols == matC%ncols)
!
        call dgemm('T', 'N', b_m, b_n, b_k, &
                   alpha, matA%m, b_lda, matB%m, b_ldb, &
                   beta, matC%m, b_ldc)
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_dgemv_N(alpha, matA, x, beta, y)
!
        implicit none
!
        type(HHO_matrix), intent(in) :: matA
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        real(kind=8), intent(in) :: alpha, beta
!
! --------------------------------------------------------------------------------------------------
!
!   DGEMM blas routine for hho_matrix
! --------------------------------------------------------------------------------------------------
!
        blas_int :: b_n, b_m, b_lda
        blas_int, parameter :: one = to_blas_int(1)
!
        b_lda = to_blas_int(matA%max_nrows)
        b_m = to_blas_int(matA%nrows)
        b_n = to_blas_int(matA%ncols)
        call dgemv('N', b_m, b_n, alpha, matA%m, &
                   b_lda, x, one, beta, y, one)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_dgemv_T(alpha, matA, x, beta, y)
!
        implicit none
!
        type(HHO_matrix), intent(in) :: matA
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        real(kind=8), intent(in) :: alpha, beta
!
! --------------------------------------------------------------------------------------------------
!
!   DGEMM blas routine for hho_matrix
! --------------------------------------------------------------------------------------------------
!
        blas_int :: b_n, b_m, b_lda
        blas_int, parameter :: one = to_blas_int(1)
!
        b_lda = to_blas_int(matA%max_nrows)
        b_m = to_blas_int(matA%nrows)
        b_n = to_blas_int(matA%ncols)
        call dgemv('T', b_m, b_n, alpha, matA%m, &
                   b_lda, x, one, beta, y, one)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hho_dsymv_U(alpha, matA, x, beta, y)
!
        implicit none
!
        type(HHO_matrix), intent(in) :: matA
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        real(kind=8), intent(in) :: alpha, beta
!
! --------------------------------------------------------------------------------------------------
!
!   DGEMM blas routine for hho_matrix
! --------------------------------------------------------------------------------------------------
!
        blas_int :: b_n, b_lda
        blas_int, parameter :: one = to_blas_int(1)
!
        b_lda = to_blas_int(matA%max_nrows)
        b_n = to_blas_int(matA%nrows)
!
        ASSERT(matA%nrows == matA%ncols)
        call dsymv('U', b_n, alpha, matA%m, &
                   b_lda, x, one, beta, y, one)
!
    end subroutine
!
end module
