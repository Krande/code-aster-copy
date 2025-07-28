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
    subroutine majou(model, modmec, solveu, num, nu,&
                     ma, mateco, moint, ndble, icor,&
                     tabad)
        character(len=2) :: model
        character(len=8) :: modmec
        character(len=19) :: solveu
        character(len=14) :: num
        character(len=14) :: nu
        character(len=8) :: ma
        character(len=*) :: mateco
        character(len=8) :: moint
        integer(kind=8) :: ndble
        integer(kind=8) :: icor(2)
        integer(kind=8) :: tabad(5)
    end subroutine majou
end interface
