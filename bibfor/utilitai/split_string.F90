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

subroutine split_string(instring, delim, string1, string2)
    ! split a string into 2 either side of a delimiter token
    !
    ! Arguments:
    !   instring (str) : input string
    !   delim (str): delimiter
    !
    ! Returns:
    !   string1 (str): first part, before the delimiter
    !   string2 (str, optional): second part, after the delimiter
    !
    implicit none
    character(len=*), intent(in) :: instring
    character(len=1), intent(in) :: delim
    character(len=*), intent(out):: string1
    character(len=*), intent(out), optional:: string2

    character(len=80) :: instr, part1, part2
    integer(kind=8) :: index

    instr = instring
    instr = trim(instr)

    index = scan(instr, delim)
    part1 = instr
    part2 = " "
    if (index .gt. 0) then
        part1 = instr(1:index-1)
        part2 = instr(index+1:)
    end if
    string1 = part1
    if (present(string2)) then
        string2 = part2
    end if
end subroutine split_string
