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
    subroutine trmult(modsta, numexi, mailla, neq, iddeeq,&
                      pside, numddl)
        integer(kind=8) :: neq
        character(len=8) :: modsta
        integer(kind=8) :: numexi
        character(len=8) :: mailla
        integer(kind=8) :: iddeeq
        real(kind=8) :: pside(neq)
        character(len=8) :: numddl
    end subroutine trmult
end interface
