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
    subroutine lcpima(fami, kpg, ksp, poum, mate,&
                      compor, instam, instap, crit, sigm,&
                      vim)
        character(len=*) :: fami
        integer(kind=8) :: kpg
        integer(kind=8) :: ksp
        character(len=1) :: poum
        integer(kind=8) :: mate
        character(len=16) :: compor
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: crit(*)
        real(kind=8) :: sigm(*)
        real(kind=8) :: vim(*)
    end subroutine lcpima
end interface
