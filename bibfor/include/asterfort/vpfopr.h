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
    subroutine vpfopr(option, typres, lmasse, lraide, ldynam,&
                      omemin, omemax, omeshi, nbfreq, npivot,&
                      omecor, precsh, nbrssa, nblagr, solveu,&
                      det, idet)
        character(len=*) :: option
        character(len=16) :: typres
        integer(kind=8) :: lmasse
        integer(kind=8) :: lraide
        integer(kind=8) :: ldynam
        real(kind=8) :: omemin
        real(kind=8) :: omemax
        real(kind=8) :: omeshi
        integer(kind=8) :: nbfreq
        integer(kind=8) :: npivot(2)
        real(kind=8) :: omecor
        real(kind=8) :: precsh
        integer(kind=8) :: nbrssa
        integer(kind=8) :: nblagr
        character(len=19) :: solveu
        real(kind=8) :: det(2)
        integer(kind=8) :: idet(2)
    end subroutine vpfopr
end interface
