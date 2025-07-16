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
    subroutine w039c2(nuzone, nuzonr, jvale, jdesc, nomgd,&
                      ifm, ifr)
        integer(kind=8) :: nuzone
        integer(kind=8) :: nuzonr
        integer(kind=8) :: jvale
        integer(kind=8) :: jdesc
        character(len=8) :: nomgd
        integer(kind=8) :: ifm
        integer(kind=8) :: ifr
    end subroutine w039c2
end interface
