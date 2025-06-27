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

function jvinfo(info)
    implicit none
    integer(kind=8) :: jvinfo
    integer(kind=8) :: info
! ----------------------------------------------------------------------
! Define the JEVEUX logging level (only used by jedet*)
!
! The level must be initialized at startup.
! The level is stored in a common for performance reasons.
!
! IN  info   Logging level
!            if info >= 0, use 'info" as the new level.
!            if info < 0, return the current level.
! ----------------------------------------------------------------------

    integer(kind=8) :: nivo
    common/jvnivo/nivo

    if (info .ge. 0) then
        nivo = info
    end if
    jvinfo = nivo
end function
