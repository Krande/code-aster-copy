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
    subroutine hujddd(carac, k, mater, ind, yf,&
                      vin, vec, mat, iret)
        character(len=6) :: carac
        integer(kind=8) :: k
        real(kind=8) :: mater(22, 2)
        integer(kind=8) :: ind(7)
        real(kind=8) :: yf(18)
        real(kind=8) :: vin(*)
        real(kind=8) :: vec(6)
        real(kind=8) :: mat(6, 6)
        integer(kind=8) :: iret
    end subroutine hujddd
end interface
