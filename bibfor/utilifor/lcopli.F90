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

subroutine lcopli(typ, mod, mater, hook)
    implicit none
!       OPERATEUR DE RIGIDITE POUR COMPORTEMENT ELASTIQUE LINEAIRE
!       IN  TYP    :  TYPE OPERATEUR
!                     'ISOTROPE'
!                     'ORTHOTRO'
!                     'ANISOTRO'
!           MOD    :  MODELISATION
!           MATER  :  COEFFICIENTS MATERIAU ELASTIQUE
!       OUT HOOK   :  OPERATEUR RIGIDITE ELASTIQUE LINEAIRE
!       ----------------------------------------------------------------
!
    integer(kind=8) :: ndt, ndi, i, j
    real(kind=8) :: un, d12, zero, deux
    parameter(d12=.5d0)
    parameter(un=1.d0)
    parameter(zero=0.d0)
    parameter(deux=2.d0)
!
    real(kind=8) :: hook(6, 6)
    real(kind=8) :: mater(*), e, nu, al, la, mu
!
    character(len=8) :: mod, typ
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
!       ----------------------------------------------------------------
!
    hook(:, :) = zero
!
    if (typ .eq. 'ISOTROPE') then
        e = mater(1)
        nu = mater(2)
        al = e*(un-nu)/(un+nu)/(un-deux*nu)
        la = nu*e/(un+nu)/(un-deux*nu)
        mu = e*d12/(un+nu)
!
! - 3D/DP/AX/CP
!
        if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'D_PLAN' .or. mod(1:6) .eq. 'C_PLAN' .or. &
            mod(1:4) .eq. 'AXIS') then
            do i = 1, ndi
                do j = 1, ndi
                    if (i .eq. j) hook(i, j) = al
                    if (i .ne. j) hook(i, j) = la
                end do
            end do
            do i = ndi+1, ndt
                do j = ndi+1, ndt
                    if (i .eq. j) hook(i, j) = deux*mu
                end do
            end do

!
! - 1D
!
        else if (mod(1:2) .eq. '1D') then
            hook(1, 1) = e
        end if
!
    else if (typ .eq. 'ORTHOTRO') then
!
        do i = 1, 6
            do j = 1, 6
                hook(i, j) = mater(6*(j-1)+i)
            end do
        end do
!
    end if
!
end subroutine
