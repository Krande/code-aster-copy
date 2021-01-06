! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1306
!
! ==================================================================================================
!
! Module for utilities for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Utilities_module
! ==================================================================================================
implicit none
! ==================================================================================================
public  :: prodBTDB, updateMatrSyme
! ==================================================================================================
private
#include "asterf_types.h"
#include "asterfort/assert.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! prodBTDB
!
! Compute product tBDB
!
! --------------------------------------------------------------------------------------------------
subroutine prodBTDB(matrD, nls, ncb, matrB, tBDB)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    integer, intent(in) :: nls, ncb
    real(kind=8), intent(in) :: matrD(nls, nls), matrB(nls, ncb)
    real(kind=8), intent(out) :: tBDB(ncb, ncb)
! - Local
    integer :: i, j, k
    real(kind=8) :: aux, DB(nls, ncb)
!   ------------------------------------------------------------------------------------------------
!
    tBDB = 0.d0
    do i = 1, nls
        do j = 1, ncb
            aux = 0.d0
            do k = 1, nls
                aux     = aux + matrD(i,k)*matrB(k,j)
                DB(i,j) = aux
            end do
        end do
    end do
!
    do i = 1, ncb
        do j = 1, ncb
            aux = 0.d0
            do k = 1, nls
                aux       = aux + matrB(k,i)*DB(k,j)
                tBDB(i,j) = aux
            end do
        end do
    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! updateMatrSyme
!
! Update symmetric matrix
!
! In  nbDof            : size of matrix
! In  matrToAdd        : matrix to add (with coefficient)
! IO  matr             : matrix to update
! In  coef             : coefficient
!
! --------------------------------------------------------------------------------------------------
subroutine updateMatrSyme(nbDof, matrToAdd, matr, coef_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    integer, intent(in)                :: nbDof
    real(kind=8), intent(in)           :: matrToAdd(nbDof, nbDof)
    real(kind=8), intent(inout)        :: matr(*)
    real(kind=8), optional, intent(in) :: coef_
! - Local
    integer :: i, j, k
    real(kind=8) :: coef
!   ------------------------------------------------------------------------------------------------
!
    coef = 1.d0
    if (present(coef_)) then
        coef = coef_
    endif
    k = 0
    do i = 1, nbDof
        do j = 1, i
            k = k + 1
            matr(k) = matr(k) + coef*matrToAdd(i, j)
        enddo
    enddo
!
!   ------------------------------------------------------------------------------------------------
end subroutine
!
end module SolidShell_Utilities_module
