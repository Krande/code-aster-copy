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
    subroutine foec2f(iuni, v, nbcoup, n1, n2,&
                      nompar, nomres)
        integer(kind=8) :: nbcoup
        integer(kind=8) :: iuni
        real(kind=8) :: v(2*nbcoup)
        integer(kind=8) :: n1
        integer(kind=8) :: n2
        character(len=*) :: nompar
        character(len=*) :: nomres
    end subroutine foec2f
end interface
