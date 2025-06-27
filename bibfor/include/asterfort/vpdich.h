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
    subroutine vpdich(lraide, lmasse, ldynam, tol, mxdich,&
                      mxfreq, nfreq, valp, ieme, det,&
                      idet, nbpas, typres, nblagr, solveu)
        integer(kind=8) :: lraide
        integer(kind=8) :: lmasse
        integer(kind=8) :: ldynam
        real(kind=8) :: tol
        integer(kind=8) :: mxdich
        integer(kind=8) :: mxfreq
        integer(kind=8) :: nfreq
        real(kind=8) :: valp(*)
        integer(kind=8) :: ieme(*)
        real(kind=8) :: det(*)
        integer(kind=8) :: idet(*)
        integer(kind=8) :: nbpas(*)
        character(len=16) :: typres
        integer(kind=8) :: nblagr
        character(len=19) :: solveu
    end subroutine vpdich
end interface
