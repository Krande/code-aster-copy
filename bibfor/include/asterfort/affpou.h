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
    subroutine affpou(tmp, tmpf, fcx, nom, isec,&
                      ivar, car, ncar, val, tab,&
                      exp, nbo, ioc, ier)
        character(len=24) :: tmp
        character(len=24) :: tmpf
        character(len=8) :: fcx
        character(len=24) :: nom
        integer(kind=8) :: isec
        integer(kind=8) :: ivar
        character(len=8) :: car(*)
        integer(kind=8) :: ncar
        real(kind=8) :: val(*)
        character(len=8) :: tab(*)
        character(len=8) :: exp(*)
        integer(kind=8) :: nbo
        character(len=6) :: ioc
        integer(kind=8) :: ier
    end subroutine affpou
end interface
