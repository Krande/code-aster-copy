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

subroutine lcopil(typ, mod, mater, kooh)
    implicit none
!     OPERATEUR DE SOUPLESSE POUR UN COMPORTEMENT ELASTIQUE LINEAIRE
!     IN  TYP    :  TYPE OPERATEUR
!                   'ISOTROPE'
!                   'ORTHOTRO'
!                   'ANISOTRO'
!         MOD    :  MODELISATION
!         MATER  :  COEFFICIENTS MATERIAU ELASTIQUE
!     OUT KOOH   :  OPERATEUR DE SOUPLESSE ELASTIQUE LINEAIRE
!     ----------------------------------------------------------------
!
    integer(kind=8) :: ndt, ndi, i, j
    real(kind=8) :: un, zero
    parameter(un=1.d0)
    parameter(zero=0.d0)
!
    real(kind=8) :: kooh(6, 6)
    real(kind=8) :: mater(*), e, nu, unpnue, unsure, mnuse
!
    character(len=8) :: mod, typ
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
!
    kooh(:, :) = zero
!
    if (typ .eq. 'ISOTROPE') then
        e = mater(1)
        nu = mater(2)
        unpnue = (un+nu)/e
        unsure = un/e
        mnuse = -nu/e
!
! - 3D/DP/AX/CP
!
        if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN' .or. mod(1:6) .eq. 'C_PLAN' .or. &
            mod(1:4) .eq. 'AXIS') then
            do i = 1, ndi
                do j = 1, ndi
                    if (i .eq. j) kooh(i, j) = unsure
                    if (i .ne. j) kooh(i, j) = mnuse
                end do
            end do
            do i = ndi+1, ndt
                do j = ndi+1, ndt
                    if (i .eq. j) kooh(i, j) = unpnue
                end do
            end do

!
! - 1D
!
        else if (mod(1:2) .eq. '1D') then
            kooh(1, 1) = unsure
        end if
!
    else if (typ .eq. 'ORTHOTRO') then
!
        do i = 1, 6
            do j = 1, 6
                kooh(i, j) = mater(36+6*(j-1)+i)
            end do
        end do

!
    end if
end subroutine
