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
    subroutine reci2d(lirela, mailla, nnoeca, noebe, nbcnx,&
                      cxma, normal, itria, xbar, iproj,&
                      excent)
        character(len=19) :: lirela
        character(len=8) :: mailla
        character(len=8) :: nnoeca
        integer(kind=8) :: noebe
        integer(kind=8) :: nbcnx
        integer(kind=8) :: cxma(*)
        real(kind=8) :: normal(*)
        integer(kind=8) :: itria
        real(kind=8) :: xbar(*)
        integer(kind=8) :: iproj
        real(kind=8) :: excent
    end subroutine reci2d
end interface
