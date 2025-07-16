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
    subroutine pj3dgb(ino2, geom2, geom1, tetr4, ndec,&
                      btdi, btvr, btnb, btlc, btco,&
                      p1, q1, r1, p2, q2,&
                      r2)
        integer(kind=8) :: ino2
        real(kind=8) :: geom2(*)
        real(kind=8) :: geom1(*)
        integer(kind=8) :: tetr4(*)
        integer(kind=8) :: ndec
        integer(kind=8) :: btdi(*)
        real(kind=8) :: btvr(*)
        integer(kind=8) :: btnb(*)
        integer(kind=8) :: btlc(*)
        integer(kind=8) :: btco(*)
        integer(kind=8) :: p1
        integer(kind=8) :: q1
        integer(kind=8) :: r1
        integer(kind=8) :: p2
        integer(kind=8) :: q2
        integer(kind=8) :: r2
    end subroutine pj3dgb
end interface
