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
    subroutine dgffibres(nboccfib, iinbgf, tousgroupesnom, tousgroupesnbf, maxmailgrp, &
                         ulnbnoeuds, ulnbmailles, nbfibres1, maxfibre1, ncarfi1)
        integer(kind=8) :: nboccfib
        integer(kind=8) :: iinbgf
        integer(kind=8) :: maxmailgrp
        integer(kind=8) :: ulnbnoeuds
        integer(kind=8) :: ulnbmailles
        integer(kind=8) :: nbfibres1
        integer(kind=8) :: tousgroupesnbf(*)
        integer(kind=8) :: maxfibre1
        integer(kind=8) :: ncarfi1
        character(len=24) :: tousgroupesnom(*)
    end subroutine dgffibres
end interface
