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
    subroutine exmali(basmod, nomint, numint, nommat, base,&
                      nblig, nbcol, ord, ii)
        character(len=8) :: basmod
        character(len=8) :: nomint
        integer(kind=8) :: numint
        character(len=24) :: nommat
        character(len=1) :: base
        integer(kind=8) :: nblig
        integer(kind=8) :: nbcol
        integer(kind=8) :: ord
        integer(kind=8) :: ii
    end subroutine exmali
end interface
