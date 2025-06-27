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
    subroutine lcicma(a, la, ca, lc, cc,&
                      xa, ya, b, lb, cb,&
                      xb, yb)
        integer(kind=8) :: cb
        integer(kind=8) :: lb
        integer(kind=8) :: ca
        integer(kind=8) :: la
        real(kind=8) :: a(la, ca)
        integer(kind=8) :: lc
        integer(kind=8) :: cc
        integer(kind=8) :: xa
        integer(kind=8) :: ya
        real(kind=8) :: b(lb, cb)
        integer(kind=8) :: xb
        integer(kind=8) :: yb
    end subroutine lcicma
end interface
