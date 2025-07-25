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
    subroutine analybar(x3d1, x3d2, x3d3, x3dp,&
                    xbar, excent, iproj, inoeu, icote)
        integer(kind=8) :: nbsom
        real(kind=8) :: x3d1(3)
        real(kind=8) :: x3d2(3)
        real(kind=8) :: x3d3(3)
        real(kind=8) :: x3dp(3)
        real(kind=8) :: xbar(3)
        real(kind=8) :: excent
        integer(kind=8) :: iproj
        integer(kind=8) :: inoeu
        integer(kind=8) :: icote
    end subroutine analybar
end interface
