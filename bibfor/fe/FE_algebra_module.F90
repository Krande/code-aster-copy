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
module FE_algebra_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "blas/dgemv.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE - Finite Element
!
! Module to add linear algebra operation
!
! --------------------------------------------------------------------------------------------------
!
    public :: dgemv_T_2xn, dgemv_T_3xn, dgemv_T_4xn, dgemv_T_6xn
    public :: dgemv_2x2, dgemv_3x3, dgemv_T_4x4, dgemv_T_6x6
    public :: daxpy_1, daxpy_1x2, daxpy_1x3, daxpy_1xm, dcopy_1
!
contains
!
! define to use hard coded loop or blas directly
! #define FE_USE_BLAS
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_T_6x6(mat, x, y, alpha)
!
        implicit none
!
        real(kind=8), intent(in) :: mat(6, *)
        real(kind=8), intent(in) :: x(*), alpha
        real(kind=8), intent(out) :: y(*)
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y = alpha * mat^T*x
!
!
! --------------------------------------------------------------------------------------------------
!
!
#ifdef FE_USE_BLAS
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n

        b_lda = to_blas_int(6)
        b_m = to_blas_int(6)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('T', b_m, b_n, alpha, mat, &
                   b_lda, x, b_incx, 0.0, y, &
                   b_incy)
#else
!
        y(1) = alpha*( &
               mat(1, 1)*x(1)+mat(2, 1)*x(2)+mat(3, 1)*x(3)+mat(4, 1)*x(4)+mat(5, 1)*x(5)+mat(6, &
               &1)*x(6) &
               )
        y(2) = alpha*( &
               mat(1, 2)*x(1)+mat(2, 2)*x(2)+mat(3, 2)*x(3)+mat(4, 2)*x(4)+mat(5, 2)*x(5)+mat(6, &
               &2)*x(6) &
               )
        y(3) = alpha*( &
               mat(1, 3)*x(1)+mat(2, 3)*x(2)+mat(3, 3)*x(3)+mat(4, 3)*x(4)+mat(5, 3)*x(5)+mat(6, &
               &3)*x(6) &
               )
        y(4) = alpha*( &
               mat(1, 4)*x(1)+mat(2, 4)*x(2)+mat(3, 4)*x(3)+mat(4, 4)*x(4)+mat(5, 4)*x(5)+mat(6, &
               &4)*x(6) &
               )
        y(5) = alpha*( &
               mat(1, 5)*x(1)+mat(2, 5)*x(2)+mat(3, 5)*x(3)+mat(4, 5)*x(4)+mat(5, 5)*x(5)+mat(6, &
               &5)*x(6) &
               )
        y(6) = alpha*( &
               mat(1, 6)*x(1)+mat(2, 6)*x(2)+mat(3, 6)*x(3)+mat(4, 6)*x(4)+mat(5, 6)*x(5)+mat(6, &
               &6)*x(6) &
               )
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_T_4x4(mat, x, y, alpha)
!
        implicit none
!
        real(kind=8), intent(in) :: mat(6, *)
        real(kind=8), intent(in) :: x(*), alpha
        real(kind=8), intent(out) :: y(*)
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y = alpha * mat^T*x
!
!
! --------------------------------------------------------------------------------------------------
!
!
#ifdef FE_USE_BLAS
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
        b_lda = to_blas_int(6)
        b_m = to_blas_int(4)
        b_n = to_blas_int(4)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('T', b_m, b_n, alpha, mat, &
                   b_lda, x, b_incx, 0.0, y, &
                   b_incy)
#else
!
        y(1) = alpha*(mat(1, 1)*x(1)+mat(2, 1)*x(2)+mat(3, 1)*x(3)+mat(4, 1)*x(4))
        y(2) = alpha*(mat(1, 2)*x(1)+mat(2, 2)*x(2)+mat(3, 2)*x(3)+mat(4, 2)*x(4))
        y(3) = alpha*(mat(1, 3)*x(1)+mat(2, 3)*x(2)+mat(3, 3)*x(3)+mat(4, 3)*x(4))
        y(4) = alpha*(mat(1, 4)*x(1)+mat(2, 4)*x(2)+mat(3, 4)*x(3)+mat(4, 4)*x(4))
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_3x3(mat, x, y, alpha)
!
        implicit none
!
        real(kind=8), intent(in) :: mat(3, *)
        real(kind=8), intent(in) :: x(*), alpha
        real(kind=8), intent(out) :: y(*)
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y = alpha * mat*x
!
!
! --------------------------------------------------------------------------------------------------
!
!
#ifdef FE_USE_BLAS
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
        b_lda = to_blas_int(3)
        b_m = to_blas_int(3)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, alpha, mat, &
                   b_lda, x, b_incx, 0.0, y, &
                   b_incy)
#else
!
        y(1) = alpha*(mat(1, 1)*x(1)+mat(1, 2)*x(2)+mat(1, 3)*x(3))
        y(2) = alpha*(mat(2, 1)*x(1)+mat(2, 2)*x(2)+mat(2, 3)*x(3))
        y(3) = alpha*(mat(3, 1)*x(1)+mat(3, 2)*x(2)+mat(3, 3)*x(3))
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_2x2(mat, x, y, alpha)
!
        implicit none
!
        real(kind=8), intent(in) :: mat(3, *)
        real(kind=8), intent(in) :: x(*), alpha
        real(kind=8), intent(out) :: y(*)
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y = alpha * mat*x
!
!
! --------------------------------------------------------------------------------------------------
!
!
#ifdef FE_USE_BLAS
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
        b_lda = to_blas_int(3)
        b_m = to_blas_int(2)
        b_n = to_blas_int(2)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, alpha, mat, &
                   b_lda, x, b_incx, 0.0, y, &
                   b_incy)
#else
!
        y(1) = alpha*(mat(1, 1)*x(1)+mat(1, 2)*x(2))
        y(2) = alpha*(mat(2, 1)*x(1)+mat(2, 2)*x(2))
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_T_6xn(mat, ncol, x, y, offset)
!
        implicit none
!
        real(kind=8), intent(in) :: mat(6, *)
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: ncol, offset
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += mat^T * x
!
! --------------------------------------------------------------------------------------------------
!
!
#ifdef FE_USE_BLAS
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
        b_lda = to_blas_int(6)
        b_m = to_blas_int(6)
        b_n = to_blas_int(ncol)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(offset)
        call dgemv('T', b_m, b_n, 1.d0, mat, &
                   b_lda, x, b_incx, 1.d0, y, &
                   b_incy)
#else
        integer(kind=8) :: icol, ind
        real(kind=8) :: tmp
!
        select case (ncol)
        case (3)
            ind = 1
            do icol = 1, 3
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (4)
            ind = 1
            do icol = 1, 4
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (5)
            ind = 1
            do icol = 1, 5
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (6)
            ind = 1
            do icol = 1, 6
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (7)
            ind = 1
            do icol = 1, 7
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (8)
            ind = 1
            do icol = 1, 8
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (9)
            ind = 1
            do icol = 1, 9
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (19)
            ind = 1
            do icol = 1, 19
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (26)
            ind = 1
            do icol = 1, 26
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case default
            ind = 1
            do icol = 1, ncol
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)+mat&
                      &(5, icol)*x(5)+mat(6, icol)*x(6)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        end select
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_T_4xn(mat, ncol, x, y, offset)
!
        implicit none
!
        real(kind=8), intent(in) :: mat(6, *)
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: ncol, offset
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += mat^T * x
!
! --------------------------------------------------------------------------------------------------
!
!
#ifdef FE_USE_BLAS
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
        b_lda = to_blas_int(4)
        b_m = to_blas_int(6)
        b_n = to_blas_int(ncol)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(offset)
        call dgemv('T', b_m, b_n, 1.d0, mat, &
                   b_lda, x, b_incx, 1.d0, y, &
                   b_incy)
#else
        integer(kind=8) :: icol, ind
        real(kind=8) :: tmp
!
        select case (ncol)
        case (2)
            ind = 1
            do icol = 1, 2
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (3)
            ind = 1
            do icol = 1, 3
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (4)
            ind = 1
            do icol = 1, 4
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (5)
            ind = 1
            do icol = 1, 5
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (6)
            ind = 1
            do icol = 1, 6
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (7)
            ind = 1
            do icol = 1, 7
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (8)
            ind = 1
            do icol = 1, 8
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case (9)
            ind = 1
            do icol = 1, 9
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        case default
            ind = 1
            do icol = 1, ncol
                tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
                y(ind) = y(ind)+tmp
                ind = ind+offset
            end do
        end select
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_T_3xn(mat, ncol, x, y)
!
        implicit none
!
        real(kind=8), intent(in) :: mat(3, *)
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: ncol
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += mat^T * x
!
! --------------------------------------------------------------------------------------------------
!
!
#ifdef FE_USE_BLAS
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
        b_lda = to_blas_int(3)
        b_m = to_blas_int(3)
        b_n = to_blas_int(ncol)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('T', b_m, b_n, 1.d0, mat, &
                   b_lda, x, b_incx, 1.d0, y, &
                   b_incy)
#else
        integer(kind=8) :: icol
!
        select case (ncol)
        case (4)
            do icol = 1, 4
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)
            end do
        case (5)
            do icol = 1, 5
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)
            end do
        case (6)
            do icol = 1, 6
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)
            end do
        case (8)
            do icol = 1, 8
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)
            end do
        case (10)
            do icol = 1, 10
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)
            end do
        case (20)
            do icol = 1, 20
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)
            end do
        case (27)
            do icol = 1, 27
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)
            end do
        case default
            do icol = 1, ncol
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)
            end do
        end select
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_T_2xn(mat, ncol, x, y)
!
        implicit none
!
        real(kind=8), intent(in) :: mat(3, *)
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: ncol
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += mat^T * x
!
! --------------------------------------------------------------------------------------------------
!
!
#ifdef FE_USE_BLAS
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n

        b_lda = to_blas_int(3)
        b_m = to_blas_int(2)
        b_n = to_blas_int(ncol)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('T', b_m, b_n, 1.d0, mat, &
                   b_lda, x, b_incx, 1.d0, y, b_incy)
#else
        integer(kind=8) :: icol
!
        select case (ncol)
        case (2)
            do icol = 1, 2
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        case (3)
            do icol = 1, 3
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        case (4)
            do icol = 1, 4
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        case (5)
            do icol = 1, 5
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        case (6)
            do icol = 1, 6
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        case (7)
            do icol = 1, 7
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        case (8)
            do icol = 1, 8
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        case (9)
            do icol = 1, 9
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        case default
            do icol = 1, ncol
                y(icol) = y(icol)+mat(1, icol)*x(1)+mat(2, icol)*x(2)
            end do
        end select
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine daxpy_1(n, alpha, x, y)
!
        implicit none
!
        real(kind=8), intent(in) :: alpha
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: n
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += y * alpha*x
!
! --------------------------------------------------------------------------------------------------
!
!
        blas_int :: b_incx, b_incy, b_n

#ifdef FE_USE_BLAS
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, alpha, x, b_incx, y, b_incy)
#else
!
        select case (n)
        case (1)
            y(1) = y(1)+alpha*x(1)
        case (2)
            y(1) = y(1)+alpha*x(1)
            y(2) = y(2)+alpha*x(2)
        case (3)
            y(1) = y(1)+alpha*x(1)
            y(2) = y(2)+alpha*x(2)
            y(3) = y(3)+alpha*x(3)
        case (4)
            y(1) = y(1)+alpha*x(1)
            y(2) = y(2)+alpha*x(2)
            y(3) = y(3)+alpha*x(3)
            y(4) = y(4)+alpha*x(4)
        case (6)
            y(1) = y(1)+alpha*x(1)
            y(2) = y(2)+alpha*x(2)
            y(3) = y(3)+alpha*x(3)
            y(4) = y(4)+alpha*x(4)
            y(5) = y(5)+alpha*x(5)
            y(6) = y(6)+alpha*x(6)
        case default
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, alpha, x, b_incx, y, b_incy)
        end select
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine daxpy_1x2(n, alpha, x, y)
!
        implicit none
!
        real(kind=8), intent(in) :: alpha
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: n
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += y * alpha*x
!
! --------------------------------------------------------------------------------------------------
!
!
        blas_int :: b_incx, b_incy, b_n

#ifdef FE_USE_BLAS
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(2)
        call daxpy(b_n, alpha, x, b_incx, y, b_incy)
#else
!
        select case (n)
        case (1)
            y(1) = y(1)+alpha*x(1)
        case (2)
            y(1) = y(1)+alpha*x(1)
            y(3) = y(3)+alpha*x(2)
        case (3)
            y(1) = y(1)+alpha*x(1)
            y(3) = y(3)+alpha*x(2)
            y(5) = y(5)+alpha*x(3)
        case (4)
            y(1) = y(1)+alpha*x(1)
            y(3) = y(3)+alpha*x(2)
            y(5) = y(5)+alpha*x(3)
            y(7) = y(7)+alpha*x(4)
        case (6)
            y(1) = y(1)+alpha*x(1)
            y(3) = y(3)+alpha*x(2)
            y(5) = y(5)+alpha*x(3)
            y(7) = y(7)+alpha*x(4)
            y(9) = y(9)+alpha*x(5)
            y(11) = y(11)+alpha*x(6)
        case default
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(2)
            call daxpy(b_n, alpha, x, b_incx, y, b_incy)
        end select
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine daxpy_1x3(n, alpha, x, y)
!
        implicit none
!
        real(kind=8), intent(in) :: alpha
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: n
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += y * alpha*x
!
! --------------------------------------------------------------------------------------------------
!
!
        blas_int :: b_incx, b_incy, b_n
!
#ifdef FE_USE_BLAS
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(3)
        call daxpy(b_n, alpha, x, b_incx, y, b_incy)
#else
!
        select case (n)
        case (1)
            y(1) = y(1)+alpha*x(1)
        case (2)
            y(1) = y(1)+alpha*x(1)
            y(4) = y(4)+alpha*x(2)
        case (3)
            y(1) = y(1)+alpha*x(1)
            y(4) = y(4)+alpha*x(2)
            y(7) = y(7)+alpha*x(3)
        case (4)
            y(1) = y(1)+alpha*x(1)
            y(4) = y(4)+alpha*x(2)
            y(7) = y(7)+alpha*x(3)
            y(10) = y(10)+alpha*x(4)
        case (6)
            y(1) = y(1)+alpha*x(1)
            y(4) = y(4)+alpha*x(2)
            y(7) = y(7)+alpha*x(3)
            y(10) = y(10)+alpha*x(4)
            y(13) = y(13)+alpha*x(5)
            y(16) = y(16)+alpha*x(6)
        case default
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(3)
            call daxpy(b_n, alpha, x, b_incx, y, b_incy)
        end select
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine daxpy_1xm(n, alpha, x, y, m)
!
        implicit none
!
        real(kind=8), intent(in) :: alpha
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: n, m
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += y * alpha*x
!
! --------------------------------------------------------------------------------------------------
!
!
        blas_int :: b_incx, b_incy, b_n
!
#ifdef FE_USE_BLAS
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(m)
        call daxpy(b_n, alpha, x, b_incx, y, b_incy)
#else
!
        select case (m)
        case (1)
            call daxpy_1(n, alpha, x, y)
        case (2)
            call daxpy_1x2(n, alpha, x, y)
        case (3)
            call daxpy_1x3(n, alpha, x, y)
        case default
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(m)
            call daxpy(b_n, alpha, x, b_incx, y, b_incy)
        end select
#endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dcopy_1(n, x, y)
!
        implicit none
!
        real(kind=8), intent(in) :: x(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: n
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dcopy product with given size
!   y = x
!
! --------------------------------------------------------------------------------------------------
!
!
        blas_int :: b_incx, b_incy, b_n

#ifdef FE_USE_BLAS
        b_n = to_blas_int(n)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, x, b_incx, y, b_incy)
#else
!
        select case (n)
        case (1)
            y(1) = x(1)
        case (2)
            y(1) = x(1)
            y(2) = x(2)
        case (3)
            y(1) = x(1)
            y(2) = x(2)
            y(3) = x(3)
        case (4)
            y(1) = x(1)
            y(2) = x(2)
            y(3) = x(3)
            y(4) = x(4)
        case (5)
            y(1) = x(1)
            y(2) = x(2)
            y(3) = x(3)
            y(4) = x(4)
        case (6)
            y(1) = x(1)
            y(2) = x(2)
            y(3) = x(3)
            y(4) = x(4)
            y(5) = x(5)
            y(6) = x(6)
        case default
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, x, b_incx, y, b_incy)
        end select
#endif
!
    end subroutine
!
end module
