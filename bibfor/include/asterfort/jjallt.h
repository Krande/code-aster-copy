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
    subroutine jjallt(lonoi, ic, gi, typei, ltypi,&
                      ci, jctab, jcdyn)
        integer(kind=8) :: lonoi
        integer(kind=8) :: ic
        character(len=*) :: gi
        character(len=*) :: typei
        integer(kind=8) :: ltypi
        character(len=*) :: ci
        integer(kind=8) :: jctab
        integer(kind=8) :: jcdyn
    end subroutine jjallt
end interface
