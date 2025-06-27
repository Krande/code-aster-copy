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
    subroutine calajt(j1, j, diag, col, n,&
                      itab, deb, tab, suiv, lt,&
                      ier)
        integer(kind=8) :: n
        integer(kind=8) :: j1
        integer(kind=8) :: j
        integer(kind=8) :: diag(0:n)
        integer(kind=8) :: col(*)
        integer(kind=8) :: itab
        integer(kind=8) :: deb(1:n)
        integer(kind=8) :: tab(*)
        integer(kind=8) :: suiv(*)
        integer(kind=8) :: lt
        integer(kind=8) :: ier
    end subroutine calajt
end interface
