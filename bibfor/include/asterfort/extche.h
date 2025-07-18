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
    subroutine extche(nchme2, nmaile, nummai, ncmp, nbm,&
                      nbc, indic, nssche, mcf, iocc,&
                      nbnac, nnoeud)
        character(len=19) :: nchme2
        character(len=8) :: nmaile(*)
        integer(kind=8) :: nummai(*)
        character(len=8) :: ncmp(*)
        integer(kind=8) :: nbm
        integer(kind=8) :: nbc
        character(len=6) :: indic
        character(len=19) :: nssche
        character(len=*) :: mcf
        integer(kind=8) :: iocc
        integer(kind=8) :: nbnac
        integer(kind=8) :: nnoeud(*)
    end subroutine extche
end interface
