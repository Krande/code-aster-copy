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
subroutine i2rdli(n, t, adr)
    implicit none
#include "asterf_types.h"
!
!
    integer(kind=8) :: n, t(*), adr
!
    aster_logical :: fini, trouve
    integer(kind=8) :: i, j
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    i = 1
    j = 0
!
    trouve = .false.
    fini = .false.
!
10  continue
    if ((.not. fini) .and. (i .lt. adr)) then
!
        if (t(i) .lt. n) then
!
            i = i+1
!
        else if (t(i) .eq. n) then
!
            trouve = .true.
            fini = .true.
!
        else
!
            fini = .true.
!
        end if
!
        goto 10
!
    end if
!
    if (.not. trouve) then
!
        do j = adr-1, i, -1
!
            t(j+1) = t(j)
!
        end do
!
        t(i) = n
        adr = adr+1
!
    end if
!
end subroutine
