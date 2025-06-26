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
interface
    subroutine nulili(nbLigr, listLigr, lili, base, gran_name, &
                      igds, mesh, nec, nlili, modeLocZ_)
        integer(kind=8), intent(in) :: nbLigr
        character(len=24), pointer :: listLigr(:)
        character(len=24), intent(in):: lili
        character(len=1), intent(in):: base
        character(len=8), intent(out) :: gran_name
        integer(kind=8), intent(out) :: igds
        character(len=8), intent(out) :: mesh
        integer(kind=8), intent(out) :: nec
        integer(kind=8), intent(out) :: nlili
        character(len=*), optional, intent(in) :: modeLocZ_
    end subroutine nulili
end interface
