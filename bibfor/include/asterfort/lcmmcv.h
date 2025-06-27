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
    subroutine lcmmcv(yd, dy, ddy, nr, itmax,&
                      toler, iter, r, rini, epstr,&
                      irteti)
        integer(kind=8) :: nr
        real(kind=8) :: yd(nr)
        real(kind=8) :: dy(nr)
        real(kind=8) :: ddy(nr)
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        integer(kind=8) :: iter
        real(kind=8) :: r(nr)
        real(kind=8) :: rini(nr)
        real(kind=8) :: epstr(6)
        integer(kind=8) :: irteti
    end subroutine lcmmcv
end interface
