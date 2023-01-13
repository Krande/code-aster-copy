! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine codhex(number, align, outstr, error)
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"

! Convert an integer to a string using the hexadecimal format.
!
!   number: Number to be converted (its absolute value is converted).
!   align: 'G' align on left, 'D' on right, 'D0' on right filled with '0'.
!   outstr: Output string.
!   error: If error is true and the output string is not long enough, a DVP_1
!       error is emitted. Otherwise the string is filled by '*'.

    integer, intent(in) :: number
    character(len=*), intent(in) :: align
    character(len=*), intent(out) :: outstr
    aster_logical, optional :: error

    integer :: value, length, need, offset
    character(len=1) :: fillwith
    character(len=6) :: buffer
    aster_logical :: failure

    failure = ASTER_FALSE
    if (present(error)) then
        failure = error
    end if
    fillwith = ' '
    if (align .eq. 'D0') then
        fillwith = '0'
    end if

    value = abs(number)
    if (value .ne. 0) then
        need = int(log(1.*value)/log(16.0))+1
    else
        need = 1
    end if
    length = len(outstr)
    if (need .gt. length) then
        ASSERT(.not. failure)
        outstr(1:length) = '*'
        goto 999
    end if

    write (buffer, "(Z6)") value
    offset = 0
    if (align .ne. 'G') then
        offset = length-need
        if (align .eq. 'C') then
            offset = offset/2
        end if
        outstr(1:offset) = fillwith
    end if
    outstr(offset+1:length) = buffer(6-need+1:6)

999 continue

end subroutine
