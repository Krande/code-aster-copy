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
    subroutine calcdq(proj, nub, nu, d, pqx,&
                      pqy, pqz, dq)
        integer(kind=8) :: proj
        real(kind=8) :: nub
        real(kind=8) :: nu
        real(kind=8) :: d(6, 6)
        real(kind=8) :: pqx(4)
        real(kind=8) :: pqy(4)
        real(kind=8) :: pqz(4)
        real(kind=8) :: dq(72)
    end subroutine calcdq
end interface
