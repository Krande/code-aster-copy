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
function char8_to_int(to_convert)
!
    implicit none
!
    character(len=8), intent(in) :: to_convert
    integer :: char8_to_int
    if (to_convert(1:1) .eq. 'M' .or. to_convert(1:1) .eq. 'N') then
        read (to_convert(2:8), *) char8_to_int
    else
        read (to_convert(1:8), *) char8_to_int
    end if
!
end function
