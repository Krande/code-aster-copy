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
    subroutine oriem0(kdim, type, coor, lino1, nbno1,&
                      lino2, nbno2, lino3, nbno3, ipos,&
                      indmai)
        character(len=2), intent(in) :: kdim
        character(len=8), intent(in) :: type
        real(kind=8), intent(in) :: coor(*)
        integer(kind=8), intent(in) :: lino1(*)
        integer(kind=8), intent(in) :: nbno1
        integer(kind=8), intent(in) :: lino2(*)
        integer(kind=8), intent(in) :: nbno2
        integer(kind=8), intent(in) :: lino3(*)
        integer(kind=8), intent(in) :: nbno3
        integer(kind=8), intent(out) :: ipos
        integer(kind=8), intent(out) :: indmai
    end subroutine oriem0
end interface
