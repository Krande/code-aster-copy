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

subroutine couplagpf3d(a, b, ngf, na, avean, &
                       nc, dgfa_ds, deltam, kmve66)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
!   complement de la matrice de couplage pour la plasticite
!   autres couplages : plasticite -> fluage

!   ************************************************************************
    implicit none

    integer(kind=8), intent(in) :: ngf, na, nc
    real(kind=8), intent(inout) :: a(ngf, ngf+1), b(ngf)
    real(kind=8), intent(in) :: avean, dgfa_ds(nc, 6), deltam, kmve66(6, 6)
!
    integer(kind=8) :: i, j, k
!
    do i = 1, 6
        do j = 1, na
!           couplage plasticite ->  kelvin
            a(i, 12+j) = avean*dgfa_ds(j, i)
!           couplage plasticite -> maxwell
!           attention kmve66 doit être dans la bonne base...
            a(i+6, 12+j) = dgfa_ds(j, i)*(deltam+kmve66(i, i))
            do k = 1, 6
                if (k .ne. i) then
                    a(i+6, 12+j) = a(i+6, 12+j)+dgfa_ds(j, k)*kmve66(i, k)
                end if
            end do
        end do
!       mise a zero des seconds membre de fluage (deja fait normalement)
        b(i) = 0.d0
        b(i+6) = 0.d0
    end do

end subroutine
