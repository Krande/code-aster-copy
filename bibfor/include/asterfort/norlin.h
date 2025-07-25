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
    subroutine norlin(typma, l, knumai, coor, dfonc,&
                      in, prec, a, b, c)
        character(len=3) :: typma
        integer(kind=8) :: l
        character(len=8) :: knumai
        real(kind=8) :: coor(3, *)
        real(kind=8) :: dfonc(*)
        integer(kind=8) :: in
        real(kind=8) :: prec
        real(kind=8) :: a
        real(kind=8) :: b
        real(kind=8) :: c
    end subroutine norlin
end interface
