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
#include "asterf_types.h"
!
interface
    subroutine pacou3(xold, fold, g, p, x,&
                      f, fvec, stpmax, check, tolx,&
                      vecr1, vecr2, typflu, vecr3, amor,&
                      masg, vecr4, vecr5, veci1, vg,&
                      indic, nbm, nmode, n)
        real(kind=8) :: xold(*)
        real(kind=8) :: fold
        real(kind=8) :: g(*)
        real(kind=8) :: p(*)
        real(kind=8) :: x(*)
        real(kind=8) :: f
        real(kind=8) :: fvec(*)
        real(kind=8) :: stpmax
        aster_logical :: check
        real(kind=8) :: tolx
        real(kind=8) :: vecr1(*)
        real(kind=8) :: vecr2(*)
        character(len=8) :: typflu
        real(kind=8) :: vecr3(*)
        real(kind=8) :: amor(*)
        real(kind=8) :: masg(*)
        real(kind=8) :: vecr4(*)
        real(kind=8) :: vecr5(*)
        integer(kind=8) :: veci1(*)
        real(kind=8) :: vg
        integer(kind=8) :: indic
        integer(kind=8) :: nbm
        integer(kind=8) :: nmode
        integer(kind=8) :: n
    end subroutine pacou3
end interface
