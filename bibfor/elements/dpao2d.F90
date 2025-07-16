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
subroutine dpao2d(angl_naut, irep, matr_tran)
!
    implicit none
!
!
    real(kind=8), intent(in) :: angl_naut(3)
    integer(kind=8), intent(out) :: irep
    real(kind=8), intent(out) :: matr_tran(4, 4)
!
! --------------------------------------------------------------------------------------------------
!
! Elasticity
!
! Construct transition matrix for orthotropic elasticity - 2D case
!
! --------------------------------------------------------------------------------------------------
!
! In  angl_naut        : nautical angles for define reference frame (AFFE_CARA_ELEM/MASSIF)
! Out irep             : 0 if matrix is trivial (identity), 1 otherwise
! Out matr_tran        : transition matrix
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: angl, cosa, sina
    real(kind=8), parameter :: zero = 0.d0
    real(kind=8), parameter :: un = 1.d0
    real(kind=8), parameter :: deux = 2.d0
!
! --------------------------------------------------------------------------------------------------
!
    irep = 0
    angl = angl_naut(1)
!
    if (angl .eq. zero) then
        irep = 0
    else
        cosa = cos(angl)
        sina = sin(angl)
        irep = 1
!
        matr_tran(1, 1) = cosa*cosa
        matr_tran(2, 1) = sina*sina
        matr_tran(3, 1) = zero
        matr_tran(4, 1) = -deux*cosa*sina
!
        matr_tran(1, 2) = sina*sina
        matr_tran(2, 2) = cosa*cosa
        matr_tran(3, 2) = zero
        matr_tran(4, 2) = deux*sina*cosa
!
        matr_tran(1, 3) = zero
        matr_tran(2, 3) = zero
        matr_tran(3, 3) = un
        matr_tran(4, 3) = zero
!
        matr_tran(1, 4) = sina*cosa
        matr_tran(2, 4) = -sina*cosa
        matr_tran(3, 4) = zero
        matr_tran(4, 4) = cosa*cosa-sina*sina
!
    end if
!
end subroutine
