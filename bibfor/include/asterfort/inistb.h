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
    subroutine inistb(maxnod, nbtyma, nomail, indic, permut,&
                      limail, indicf, permuf, maxfa)
        integer(kind=8) :: maxfa
        integer(kind=8) :: maxnod
        integer(kind=8) :: nbtyma
        character(len=8) :: nomail(*)
        integer(kind=8) :: indic(*)
        integer(kind=8) :: permut(maxnod, *)
        integer(kind=8) :: limail(*)
        integer(kind=8) :: indicf(*)
        integer(kind=8) :: permuf(maxfa, *)
    end subroutine inistb
end interface
