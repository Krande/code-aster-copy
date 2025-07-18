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
subroutine dmdepv(rho, satur, tbiot, dmdeps)
!
    implicit none
!
    real(kind=8), intent(in) :: rho, tbiot(6), satur
    real(kind=8), intent(out) :: dmdeps(6)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Derivative of quantity of mass by volumic mass - Mechanical part (strains)
!
! --------------------------------------------------------------------------------------------------
!
! In  rho              : volumic mass
! In  tbiot            : tensor of Biot
! In  satur            : value of saturation
! Out dmdeps           : derivative of quantity of mass by volumic mass - Mechanical part (strains)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
!
! --------------------------------------------------------------------------------------------------
!
    dmdeps(:) = 0.d0
    do i = 1, 3
        dmdeps(i) = rho*tbiot(i)*satur
    end do
    do i = 4, 6
        dmdeps(i) = rho*tbiot(i)*satur*rac2
    end do
!
end subroutine
