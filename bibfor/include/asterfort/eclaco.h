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
    subroutine eclaco(ipg, mxnbn2, connx, nbno2, i1,&
                      i2, i3, i4, i5, i6,&
                      i7, i8)
        integer(kind=8) :: mxnbn2
        integer(kind=8) :: ipg
        integer(kind=8) :: connx(mxnbn2, *)
        integer(kind=8) :: nbno2(*)
        integer(kind=8) :: i1
        integer(kind=8) :: i2
        integer(kind=8) :: i3
        integer(kind=8) :: i4
        integer(kind=8) :: i5
        integer(kind=8) :: i6
        integer(kind=8) :: i7
        integer(kind=8) :: i8
    end subroutine eclaco
end interface
