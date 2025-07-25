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
    subroutine irgme2(numold, ima, connex, nbord2, tabd,&
                      tabl, tabv, partie, jtype, nbno,&
                      listno, nbcmp, ifi, iadmax)
        integer(kind=8) :: numold(*)
        integer(kind=8) :: ima
        character(len=24) :: connex
        integer(kind=8) :: nbord2
        integer(kind=8) :: tabd(*)
        integer(kind=8) :: tabl(*)
        integer(kind=8) :: tabv(*)
        character(len=*) :: partie
        integer(kind=8) :: jtype
        integer(kind=8) :: nbno
        integer(kind=8) :: listno(*)
        integer(kind=8) :: nbcmp
        integer(kind=8) :: ifi
        integer(kind=8) :: iadmax
    end subroutine irgme2
end interface
