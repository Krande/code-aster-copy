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
    subroutine statpu(nbobst, nbpt, temps, fcho, vgli,&
                      iadh, wk1, wk2, wk3, iwk4,&
                      idebut, nbloc, nbval, ifires, inoe,&
                      impr, pusurn)
        integer(kind=8) :: nbobst
        integer(kind=8) :: nbpt
        real(kind=8) :: temps(*)
        real(kind=8) :: fcho(*)
        real(kind=8) :: vgli(*)
        integer(kind=8) :: iadh(*)
        real(kind=8) :: wk1(*)
        real(kind=8) :: wk2(*)
        real(kind=8) :: wk3(*)
        integer(kind=8) :: iwk4(*)
        integer(kind=8) :: idebut
        integer(kind=8) :: nbloc
        integer(kind=8) :: nbval
        integer(kind=8) :: ifires
        integer(kind=8) :: inoe
        integer(kind=8) :: impr
        real(kind=8) :: pusurn
    end subroutine statpu
end interface
