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
    subroutine exprli(basmdz, lintfz, nmintz, numint, famprz,&
                      ii, ordo)
        character(len=*) :: basmdz
        character(len=*) :: lintfz
        character(len=*) :: nmintz
        integer(kind=8) :: numint
        character(len=*) :: famprz
        integer(kind=8) :: ii
        integer(kind=8) :: ordo
    end subroutine exprli
end interface
