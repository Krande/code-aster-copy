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
    subroutine sspace(lraid, lmatra, lmass, neq, nbvec,&
                      nfreq, lprod, itemax, nperm, tol,&
                      toldyn, vect, valpro, nitjac, nitbat,&
                      solveu)
        integer(kind=8) :: nbvec
        integer(kind=8) :: neq
        integer(kind=8) :: lraid
        integer(kind=8) :: lmatra
        integer(kind=8) :: lmass
        integer(kind=8) :: nfreq
        integer(kind=8) :: lprod(neq)
        integer(kind=8) :: itemax
        integer(kind=8) :: nperm
        real(kind=8) :: tol
        real(kind=8) :: toldyn
        real(kind=8) :: vect(neq, nbvec)
        real(kind=8) :: valpro(nbvec)
        integer(kind=8) :: nitjac
        integer(kind=8) :: nitbat
        character(len=19) :: solveu
    end subroutine sspace
end interface
