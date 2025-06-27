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
    subroutine caldis(fremax, fremin, pas, frexci, nbptmd,&
                      nbmode, lismod, fremod, amomod, nindex,&
                      npdsc3, frefin)
        real(kind=8) :: fremax
        real(kind=8) :: fremin
        real(kind=8) :: pas
        character(len=4) :: frexci
        integer(kind=8) :: nbptmd
        integer(kind=8) :: nbmode
        integer(kind=8) :: lismod(*)
        real(kind=8) :: fremod(*)
        real(kind=8) :: amomod(*)
        integer(kind=8) :: nindex
        integer(kind=8) :: npdsc3
        real(kind=8) :: frefin
    end subroutine caldis
end interface
