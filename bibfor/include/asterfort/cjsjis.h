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
    subroutine cjsjis(mod, mater, deps, yd, yf,&
                      r, drdy)
        integer(kind=8), parameter :: nmod=8
        character(len=8) :: mod
        real(kind=8) :: mater(14, 2)
        real(kind=8) :: deps(6)
        real(kind=8) :: yd(nmod)
        real(kind=8) :: yf(nmod)
        real(kind=8) :: r(nmod)
        real(kind=8) :: drdy(nmod, nmod)
    end subroutine cjsjis
end interface
