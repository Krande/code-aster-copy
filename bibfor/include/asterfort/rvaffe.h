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
    subroutine rvaffe(mcf, iocc, sdlieu, sdeval, sdmail,&
                      typaff, quant, option, rep, nomtab,&
                      xnovar, ncheff, i1, isd)
        character(len=*) :: mcf
        integer(kind=8) :: iocc
        character(len=24) :: sdlieu
        character(len=19) :: sdeval
        character(len=24) :: sdmail
        character(len=1) :: typaff
        character(len=*) :: quant
        character(len=*) :: option
        character(len=*) :: rep
        character(len=19) :: nomtab
        character(len=24) :: xnovar
        character(len=16) :: ncheff
        integer(kind=8) :: i1
        integer(kind=8) :: isd
    end subroutine rvaffe
end interface
