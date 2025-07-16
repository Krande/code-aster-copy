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
subroutine moyte2(fami, npg, poum, temp, iret)
    implicit none
#include "asterfort/rcvarc.h"
    integer(kind=8) :: npg, iret, iretm
    real(kind=8) :: temp
    character(len=*) :: fami, poum
    integer(kind=8) :: kpg
    real(kind=8) :: tg, tgtot, tm, ti, ts
!
    tgtot = 0.d0
    do kpg = 1, npg
        call rcvarc(' ', 'TEMP', poum, fami, kpg, &
                    1, tm, iretm)
        call rcvarc(' ', 'TEMP_INF', poum, fami, kpg, &
                    1, ti, iret)
        call rcvarc(' ', 'TEMP_SUP', poum, fami, kpg, &
                    1, ts, iret)
        if (iret .eq. 1) then
            temp = tm
            goto 999
        else if (iretm .ne. 0) then
            tm = (ti+ts)/2.d0
        end if
        tg = (4.d0*tm+ti+ts)/6.0d0
        tgtot = tgtot+tg
    end do
    temp = tgtot/npg
999 continue
end subroutine
