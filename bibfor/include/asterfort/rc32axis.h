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
    subroutine rc32axis(nbabsc, absc, xcoo, ycoo, vale, momen0_axis, momen1_axis, momen2_axis,rho)
        integer(kind=8) :: nbabsc
        real(kind=8) :: absc(nbabsc)
        real(kind=8) :: vale(4,nbabsc)
        real(kind=8) :: xcoo(nbabsc)
        real(kind=8) :: ycoo(nbabsc)
        real(kind=8) :: momen0_axis(4)
        real(kind=8) :: momen1_axis(4)
        real(kind=8) :: momen2_axis(4)
        real(kind=8) :: rho
    end subroutine rc32axis
end interface
!
