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
    subroutine xcedge(ndime, pinref, pi1, pi2, pmiref, m12, crit)
        integer(kind=8) :: ndime
        integer(kind=8) :: pi1
        integer(kind=8) :: pi2
        integer(kind=8) :: m12
        real(kind=8) :: pinref(*)
        real(kind=8) :: pmiref(*)
        real(kind=8) :: crit
    end subroutine xcedge
end interface 
