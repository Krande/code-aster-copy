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
    public :: dgemv_T_4xn, dgemv_T_6xn
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
!
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
!
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
