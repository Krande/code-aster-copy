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
interface
    subroutine vff2dn(ndim, nno, ipg, ipoids, idfde,&
                      coor, nx, ny, jac)
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: nno
        integer(kind=8), intent(in) :: ipg
        integer(kind=8), intent(in) :: ipoids
        integer(kind=8), intent(in) :: idfde
        real(kind=8), intent(in) :: coor(*)
        real(kind=8), intent(out) :: nx
        real(kind=8), intent(out) :: ny
        real(kind=8), intent(out) :: jac
    end subroutine vff2dn
end interface
