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
    subroutine glrcmm(zimat, matr, ep, surfgp, p,&
                      epst, deps, dsig, ecr, delas,&
                      dsidep, crit, codret)
        integer(kind=8) :: zimat
        real(kind=8) :: matr(*)
        real(kind=8) :: ep
        real(kind=8) :: surfgp
        real(kind=8) :: p(3, 3)
        real(kind=8) :: epst(*)
        real(kind=8) :: deps(*)
        real(kind=8) :: dsig(*)
        real(kind=8) :: ecr(*)
        real(kind=8) :: delas(6, *)
        real(kind=8) :: dsidep(6, *)
        real(kind=8) :: crit(*)
        integer(kind=8) :: codret
    end subroutine glrcmm
end interface
