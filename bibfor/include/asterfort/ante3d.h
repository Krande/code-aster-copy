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
    subroutine ante3d(nbsom, itetra, xbar, ksi1, ksi2,&
                      ksi3)
        integer(kind=8) :: nbsom
        integer(kind=8) :: itetra
        real(kind=8) :: xbar(*)
        real(kind=8) :: ksi1
        real(kind=8) :: ksi2
        real(kind=8) :: ksi3
    end subroutine ante3d
end interface
