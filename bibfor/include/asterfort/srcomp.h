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
    subroutine srcomp(mod, imate, instam, instap,&
                      deps, sigm, vinm,&
                      option, sigp, vinp, dside, retcom,&
                      invi)
        integer(kind=8) :: invi
        character(len=8) :: mod(*)
        integer(kind=8) :: imate
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: deps(6)
        real(kind=8) :: sigm(6)
        real(kind=8) :: vinm(invi)
        character(len=16) :: option
        real(kind=8) :: sigp(6)
        real(kind=8) :: vinp(invi)
        real(kind=8) :: dside(6, 6)
        integer(kind=8) :: retcom
    end subroutine srcomp
end interface
