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

subroutine couplagfp3d(a, ngf, na, nc, &
                       dpfa_ds, dpfa_dpg, dpg_depsa6, raideur66p)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
!   autres couplages : fluage -> plasticite

!   ************************************************************************
    implicit none

    integer(kind=8), intent(in) :: ngf, na, nc
    real(kind=8), intent(inout) :: a(ngf, ngf+1)
    real(kind=8), intent(in) :: dpfa_ds(nc, 6), dpfa_dpg(nc)
    real(kind=8), intent(in) :: dpg_depsa6(6), raideur66p(6, 6)
!
    integer(kind=8) :: i, j, k
!
    do i = 1, na
        do j = 1, 6
!           couplage kelvin -> plasticite
            a(12+i, j) = 0.d0
            do k = 1, 6
                a(12+i, j) = a(12+i, j)-dpfa_ds(i, k)*raideur66p(k, j)
            end do
!           ici les deformation anelastique sont du fluage donc
!           on prend dpg_depsa6 du cas general
            a(12+i, j) = a(12+i, j)+dpfa_dpg(i)*dpg_depsa6(j)
!           couplage maxwell -> plasticite
            a(12+i, j+6) = a(12+i, j)
        end do
    end do

end subroutine
