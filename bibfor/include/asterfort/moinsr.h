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
    subroutine moinsr(j, n, idil, idiich, idsuiv,&
                      nosuiv, idip, noip, iilib, iimax)
        integer(kind=8) :: j
        integer(kind=8) :: n
        integer(kind=8) :: idil
        integer(kind=8) :: idiich
        integer(kind=8) :: idsuiv
        character(len=*) :: nosuiv
        integer(kind=8) :: idip
        character(len=*) :: noip
        integer(kind=8) :: iilib
        integer(kind=8) :: iimax
    end subroutine moinsr
end interface
