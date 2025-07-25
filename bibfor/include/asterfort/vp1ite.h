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
    subroutine vp1ite(lmasse, lraide, ldynam, x, imode,&
                      valp, neq, mxiter, tol, iter,&
                      x0, mx, err, iexcl, place,&
                      iquoti, solveu)
        integer(kind=8) :: neq
        integer(kind=8) :: lmasse
        integer(kind=8) :: lraide
        integer(kind=8) :: ldynam
        real(kind=8) :: x(neq, 1)
        integer(kind=8) :: imode
        real(kind=8) :: valp
        integer(kind=8) :: mxiter
        real(kind=8) :: tol
        integer(kind=8) :: iter
        real(kind=8) :: x0(neq)
        real(kind=8) :: mx(neq, *)
        real(kind=8) :: err
        integer(kind=8) :: iexcl(*)
        integer(kind=8) :: place
        integer(kind=8) :: iquoti
        character(len=19) :: solveu
    end subroutine vp1ite
end interface
