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
    subroutine vp2ini(ldynam, lmasse, ldynfa, neq, nbvect,&
                      nborto, prorto, ddlexc, ddllag, alpha,&
                      beta, signes, vect, prsudg, nstoc,&
                      omeshi, solveu)
        integer(kind=8) :: neq
        integer(kind=8) :: ldynam
        integer(kind=8) :: lmasse
        integer(kind=8) :: ldynfa
        integer(kind=8) :: nbvect
        integer(kind=8) :: nborto
        real(kind=8) :: prorto
        integer(kind=8) :: ddlexc(*)
        integer(kind=8) :: ddllag(*)
        real(kind=8) :: alpha(*)
        real(kind=8) :: beta(*)
        real(kind=8) :: signes(*)
        real(kind=8) :: vect(neq, *)
        real(kind=8) :: prsudg
        integer(kind=8) :: nstoc
        real(kind=8) :: omeshi
        character(len=19) :: solveu
    end subroutine vp2ini
end interface
