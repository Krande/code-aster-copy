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

interface
    subroutine mgis_integrate(extern_addr, stress, statev, ddsdde, dtime, rdt, &
                              pnewdt, retcode)
        character(len=16), intent(in) :: extern_addr
        real(kind=8), intent(inout) :: stress(*)
        real(kind=8), intent(inout) :: statev(*)
        real(kind=8), intent(out) :: ddsdde(*)
        real(kind=8), intent(in) :: rdt, dtime
        real(kind=8), intent(out) :: pnewdt
        integer(kind=8), intent(out) :: retcode
    end subroutine mgis_integrate
end interface
