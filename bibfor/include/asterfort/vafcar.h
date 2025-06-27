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
    subroutine vafcar(tpgz, imclf, nmobjz, nutyel, ntyele, car, ncar, ivr, kioc, ier)
        character(len=*) :: tpgz
        integer(kind=8) :: imclf
        character(len=*) :: nmobjz
        integer(kind=8) :: nutyel
        integer(kind=8) :: ntyele(*)
        character(len=*) :: car(*)
        integer(kind=8) :: ncar
        integer(kind=8) :: ivr(*)
        character(len=6) :: kioc
        integer(kind=8) :: ier
    end subroutine vafcar
end interface
