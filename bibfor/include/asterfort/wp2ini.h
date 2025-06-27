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
    subroutine wp2ini(appr, lmasse, lamor, lraide, lmatra,&
                      lmtpsc, sigma, xh, xb, optiof,&
                      prorto, nborto, nbvect, neq, lbloq,&
                      lddl, alpha, beta, signe, yh,&
                      yb, solveu)
        integer(kind=8) :: neq
        character(len=1) :: appr
        integer(kind=8) :: lmasse
        integer(kind=8) :: lamor
        integer(kind=8) :: lraide
        integer(kind=8) :: lmatra
        integer(kind=8) :: lmtpsc
        complex(kind=8) :: sigma
        real(kind=8) :: xh(*)
        real(kind=8) :: xb(*)
        character(len=*) :: optiof
        real(kind=8) :: prorto
        integer(kind=8) :: nborto
        integer(kind=8) :: nbvect
        integer(kind=8) :: lbloq(*)
        integer(kind=8) :: lddl(*)
        real(kind=8) :: alpha(*)
        real(kind=8) :: beta(*)
        real(kind=8) :: signe(*)
        real(kind=8) :: yh(neq, *)
        real(kind=8) :: yb(neq, *)
        character(len=19) :: solveu
    end subroutine wp2ini
end interface
