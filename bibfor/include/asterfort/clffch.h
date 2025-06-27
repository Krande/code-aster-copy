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
    subroutine clffch(alias, type, nno, xi, yi,&
                      zi, xin, yin, zin, tn,&
                      ajx, ajy, ajz, bjxx, bjyy,&
                      bjzz, bjxy, bjxz, bjyz, ider)
        character(len=6) :: alias
        character(len=6) :: type
        integer(kind=8) :: nno
        real(kind=8) :: xi
        real(kind=8) :: yi
        real(kind=8) :: zi
        real(kind=8) :: xin(20)
        real(kind=8) :: yin(20)
        real(kind=8) :: zin(20)
        real(kind=8) :: tn(*)
        real(kind=8) :: ajx(*)
        real(kind=8) :: ajy(*)
        real(kind=8) :: ajz(*)
        real(kind=8) :: bjxx(*)
        real(kind=8) :: bjyy(*)
        real(kind=8) :: bjzz(*)
        real(kind=8) :: bjxy(*)
        real(kind=8) :: bjxz(*)
        real(kind=8) :: bjyz(*)
        integer(kind=8) :: ider
    end subroutine clffch
end interface
