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
module HHO_matrix_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/readMatrix.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - dynamic matrix
!
! --------------------------------------------------------------------------------------------------
!

    type HHO_matrix
        integer :: nrows = 0, ncols = 0
        integer :: max_nrows = 0, max_ncols = 0
        aster_logical :: is_allocated = ASTER_FALSE
! ----- array
        real(kind=8), dimension(:, :), pointer :: m
!
! ----- member function
    contains
        procedure, pass :: initialize => hhoMatriceInit
        procedure, pass :: free => hhoMatriceFree
        procedure, pass :: write => hhoMatriceWrite
        procedure, pass :: read => hhoMatriceRead
        procedure, pass :: setValue => hhoMatriceSetValue
        procedure, pass :: print => hhoMatricePrint
        procedure, pass :: copySymU => hhoMatriceCopySymU
!
    end type HHO_matrix
!
    public   :: HHO_matrix
    private  :: hhoMatriceInit, hhoMatriceFree, hhoMatriceWrite, hhoMatriceSetValue
    private  :: hhoMatriceRead, hhoMatricePrint, hhoMatriceCopySymU
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
    subroutine hhoMatriceInit(this, n_rows, n_cols, val)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        integer, intent(in) :: n_rows, n_cols
        real(kind=8), intent(in), optional :: val
!
        this%nrows = n_rows
        this%ncols = n_cols
        this%max_nrows = n_rows
        this%max_ncols = n_cols
!
        ASSERT(.not. this%is_allocated)
        ASSERT(n_rows > 0 .and. n_cols > 0)
!
        allocate (this%m(n_rows, n_cols))
        this%is_allocated = ASTER_TRUE
!
        if (present(val)) then
            call this%setValue(val)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceFree(this)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
!
        this%nrows = 0
        this%ncols = 0

        if (this%is_allocated) then
            deallocate (this%m)
            this%is_allocated = ASTER_FALSE
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceWrite(this, name, l_sym)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
        character(len=*), intent(in) :: name
        aster_logical, intent(in) :: l_sym
!
        ASSERT(this%is_allocated)
        call writeMatrix(name, this%nrows, this%ncols, l_sym, this%m)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceRead(this, name, l_sym)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
        character(len=*), intent(in) :: name
        aster_logical, intent(in) :: l_sym
!
        ASSERT(this%is_allocated)
        call readMatrix(name, this%nrows, this%ncols, l_sym, this%m)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceSetValue(this, val)
!
        implicit none
!
        class(HHO_matrix), intent(inout) :: this
        real(kind=8), intent(in) :: val
!
        ASSERT(this%is_allocated)
        this%m = val
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatricePrint(this)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer :: i
!
!
        write (6, *) "matrix of", this%nrows, "rows x ", this%ncols, "cols"
        do i = 1, this%nrows
            write (6, '(50ES13.6)') this%m(i, 1:this%ncols)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMatriceCopySymU(this)
!
        implicit none
!
        class(HHO_matrix), intent(in) :: this
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer :: i
!
!
        do i = 1, this%ncols-1
            this%m(i+1:this%ncols, i) = this%m(i, i+1:this%ncols)
        end do
!
    end subroutine
!
end module
