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
    subroutine newdep(neq, c, dt, d0, v0,&
                      a0, d1, a1)
        integer(kind=8) :: neq
        real(kind=8) :: c
        real(kind=8) :: dt
        real(kind=8) :: d0(*)
        real(kind=8) :: v0(*)
        real(kind=8) :: a0(*)
        real(kind=8) :: d1(*)
        real(kind=8) :: a1(*)
    end subroutine newdep
end interface
