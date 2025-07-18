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

subroutine matrHooke3d(elas_type, angl, &
                       h, g, g1, g2, g3, &
                       matr_elas)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dpassa.h"
#include "asterfort/utbtab.h"
!
!
    integer(kind=8), intent(in) :: elas_type
    real(kind=8), intent(in) :: angl(3)
    real(kind=8), intent(in) :: g, h(6)
    real(kind=8), intent(in) :: g1, g2, g3
    real(kind=8), intent(out) :: matr_elas(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
! Hooke matrix for iso-parametric elements
!
! 3D and Fourier
!
! --------------------------------------------------------------------------------------------------
!
! In  elas_type          : Type of elasticity
!                 1 - Isotropic
!                 2 - Orthotropic
!                 3 - Transverse isotropic
!                           or viscoelasticity
!                 4 - Isotropic
!                 5 - Orthotropic
!                 6 - Transverse isotropic
! In  angl             : nautical angles
! In  h                : Hook coefficient (all)
! In  g                : shear ratio (isotropic/Transverse isotropic)
! In  g1               : shear ratio (Orthotropic)
! In  g2               : shear ratio (Orthotropic)
! In  g3               : shear ratio (Orthotropic)
! Out matr_elas        : Hooke matrix in global reference frame
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: irep, i, j
    real(kind=8) :: matr_tran(6, 6), dorth(6, 6), work(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
    matr_elas(:, :) = 0.d0
    dorth(:, :) = 0.d0
    work(:, :) = 0.d0
!
! - Compute Hooke matrix
!
    if (elas_type .eq. 1 .or. elas_type .eq. 4) then
!
! ----- Isotropic elastic matrix
!
        matr_elas(1, 1) = h(1)
        matr_elas(1, 2) = h(2)
        matr_elas(1, 3) = h(2)
        matr_elas(2, 1) = h(2)
        matr_elas(2, 2) = h(1)
        matr_elas(2, 3) = h(2)
        matr_elas(3, 1) = h(2)
        matr_elas(3, 2) = h(2)
        matr_elas(3, 3) = h(1)
        matr_elas(4, 4) = g
        matr_elas(5, 5) = g
        matr_elas(6, 6) = g
!
    else if (elas_type .eq. 2 .or. elas_type .eq. 5) then
!
! ----- Orthotropic matrix
!
        dorth(1, 1) = h(1)
        dorth(1, 2) = h(2)
        dorth(1, 3) = h(3)
        dorth(2, 2) = h(4)
        dorth(2, 3) = h(5)
        dorth(3, 3) = h(6)
        dorth(2, 1) = dorth(1, 2)
        dorth(3, 1) = dorth(1, 3)
        dorth(3, 2) = dorth(2, 3)
        dorth(4, 4) = g1
        dorth(5, 5) = g2
        dorth(6, 6) = g3
!
! ----- Compute transition matrix from orthotropic basis to global 3D basis
!
        call dpassa(angl, irep, matr_tran)
!
! ----- Change Hooke matrix to global 3D basis
!
        ASSERT((irep .eq. 1) .or. (irep .eq. 0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 6, 6, dorth, matr_tran, work, matr_elas)
        else if (irep .eq. 0) then
            do i = 1, 6
                do j = 1, 6
                    matr_elas(i, j) = dorth(i, j)
                end do
            end do
        end if
!
    else if (elas_type .eq. 3 .or. elas_type .eq. 6) then
!
! ----- Transverse isotropic matrix
!
        dorth(1, 1) = h(1)
        dorth(1, 2) = h(2)
        dorth(1, 3) = h(3)
        dorth(2, 1) = dorth(1, 2)
        dorth(2, 2) = dorth(1, 1)
        dorth(2, 3) = dorth(1, 3)
        dorth(3, 1) = dorth(1, 3)
        dorth(3, 2) = dorth(2, 3)
        dorth(3, 3) = h(4)
        dorth(4, 4) = h(5)
        dorth(5, 5) = g
        dorth(6, 6) = dorth(5, 5)
!
! ----- Compute transition matrix from orthotropic basis to global 3D basis
!
        call dpassa(angl, irep, matr_tran)
!
! ----- Change Hooke matrix to global 3D basis
!
        ASSERT((irep .eq. 1) .or. (irep .eq. 0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 6, 6, dorth, matr_tran, work, matr_elas)
        else if (irep .eq. 0) then
            do i = 1, 6
                do j = 1, 6
                    matr_elas(i, j) = dorth(i, j)
                end do
            end do
        end if
!
    else
        ASSERT(.false.)
    end if
!
end subroutine
