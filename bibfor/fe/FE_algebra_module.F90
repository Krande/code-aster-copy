! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
!
! --------------------------------------------------------------------------------------------------
!
! FE - Finite Element
!
! Module to add linear algebra operation
!
! --------------------------------------------------------------------------------------------------
!
    public :: dgemv_T_4xn, dgemv_T_6xn, dgemv_T_4x4, dgemv_T_6x6
!    private  ::
!
contains
!
! define to use hard coded loop or blas directly
#define FE_USE_BLAS 0
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine dgemv_T_6x6(mat, x, y, alpha)
!
        implicit none
!
        real(kind=8), intent(in)     :: mat(6, *)
        real(kind=8), intent(in)     :: x(*), alpha
        real(kind=8), intent(out)    :: y(*)
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
#if FE_USE_BLAS
        call dgemv('T', 6, 6, alpha, mat, 6, x, 1, 0.0, y, 1)
#else
!
        y(1) = alpha*(mat(1, 1)*x(1)+mat(2, 1)*x(2)+mat(3, 1)*x(3)+ &
                      mat(4, 1)*x(4)+mat(5, 1)*x(5)+mat(6, 1)*x(6))
        y(2) = alpha*(mat(1, 2)*x(1)+mat(2, 2)*x(2)+mat(3, 2)*x(3)+ &
                      mat(4, 2)*x(4)+mat(5, 2)*x(5)+mat(6, 2)*x(6))
        y(3) = alpha*(mat(1, 3)*x(1)+mat(2, 3)*x(2)+mat(3, 3)*x(3)+ &
                      mat(4, 3)*x(4)+mat(5, 3)*x(5)+mat(6, 3)*x(6))
        y(4) = alpha*(mat(1, 4)*x(1)+mat(2, 4)*x(2)+mat(3, 4)*x(3)+ &
                      mat(4, 4)*x(4)+mat(5, 4)*x(5)+mat(6, 4)*x(6))
        y(5) = alpha*(mat(1, 5)*x(1)+mat(2, 5)*x(2)+mat(3, 5)*x(3)+ &
                      mat(4, 5)*x(4)+mat(5, 5)*x(5)+mat(6, 5)*x(6))
        y(6) = alpha*(mat(1, 6)*x(1)+mat(2, 6)*x(2)+mat(3, 6)*x(3)+ &
                      mat(4, 6)*x(4)+mat(5, 6)*x(5)+mat(6, 6)*x(6))
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
        real(kind=8), intent(in)     :: mat(6, *)
        real(kind=8), intent(in)     :: x(*), alpha
        real(kind=8), intent(out)    :: y(*)
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
#if FE_USE_BLAS
        call dgemv('T', 4, 4, alpha, mat, 6, x, 1, 0.0, y, 1)
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
    subroutine dgemv_T_6xn(mat, ncol, x, y, offset)
!
        implicit none
!
        real(kind=8), intent(in)     :: mat(6, *)
        real(kind=8), intent(in)     :: x(*)
        real(kind=8), intent(inout)  :: y(*)
        integer, intent(in)          :: ncol, offset
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += mat * x
!
! --------------------------------------------------------------------------------------------------
!
!
#if FE_USE_BLAS
        call dgemv('T', 6, ncol, 1.d0, mat, 6, x, 1, 1.d0, y, offset)
#else
        integer :: icol, ind
        real(kind=8) :: tmp

        ind = 1
        do icol = 1, ncol
            tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+ &
                  mat(4, icol)*x(4)+mat(5, icol)*x(5)+mat(6, icol)*x(6)
            y(ind) = y(ind)+tmp
            ind = ind+offset
        end do
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
        real(kind=8), intent(in)     :: mat(6, *)
        real(kind=8), intent(in)     :: x(*)
        real(kind=8), intent(inout)  :: y(*)
        integer, intent(in)          :: ncol, offset
!
! --------------------------------------------------------------------------------------------------
!
!   Encapsulation of dgemv product with given size
!   y += mat * x
!
! --------------------------------------------------------------------------------------------------
!
!
#if FE_USE_BLAS
        call dgemv('T', 6, ncol, 1.d0, mat, 4, x, 1, 1.d0, y, offset)
#else
        integer :: icol, ind
        real(kind=8) :: tmp

        ind = 1
        do icol = 1, ncol
            tmp = mat(1, icol)*x(1)+mat(2, icol)*x(2)+mat(3, icol)*x(3)+mat(4, icol)*x(4)
            y(ind) = y(ind)+tmp
            ind = ind+offset
        end do
#endif
!
    end subroutine
!
end module
