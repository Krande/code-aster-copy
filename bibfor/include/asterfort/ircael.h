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
    subroutine ircael(jcesdi, jcesli, jcesvi, jcesci, nummai,&
                      nbqcou, nbtcou, nbrsec, nbrfib, nbrgrf, nugrfi)
        integer(kind=8) :: jcesdi
        integer(kind=8) :: jcesli
        integer(kind=8) :: jcesvi
        integer(kind=8) :: jcesci
        integer(kind=8) :: nummai
        integer(kind=8) :: nbqcou
        integer(kind=8) :: nbtcou
        integer(kind=8) :: nbrsec
        integer(kind=8) :: nbrfib
        integer(kind=8) :: nbrgrf
        integer(kind=8) :: nugrfi(10)
    end subroutine ircael
end interface
