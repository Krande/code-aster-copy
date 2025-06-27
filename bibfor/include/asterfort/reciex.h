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
    subroutine reciex(intexc, iderex, nindex, nnoeex, ncmpex,&
                      nvasex, graexc, excmod, napexc)
        character(len=8) :: intexc
        integer(kind=8) :: iderex
        integer(kind=8) :: nindex
        integer(kind=8) :: nnoeex
        integer(kind=8) :: ncmpex
        integer(kind=8) :: nvasex
        character(len=16) :: graexc
        character(len=4) :: excmod
        integer(kind=8) :: napexc
    end subroutine reciex
end interface
