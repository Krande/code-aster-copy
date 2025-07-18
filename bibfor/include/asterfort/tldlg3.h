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
    subroutine tldlg3(metrez, renum, istop, lmat, ildeb,&
                      ilfin, ndigit, ndeci, isingu, npvneg,&
                      iret, solvop)
        character(len=*) :: metrez
        character(len=*) :: renum
        integer(kind=8) :: istop
        integer(kind=8) :: lmat
        integer(kind=8) :: ildeb
        integer(kind=8) :: ilfin
        integer(kind=8) :: ndigit
        integer(kind=8) :: ndeci
        integer(kind=8) :: isingu
        integer(kind=8) :: npvneg
        integer(kind=8) :: iret
        character(len=*) :: solvop
    end subroutine tldlg3
end interface
