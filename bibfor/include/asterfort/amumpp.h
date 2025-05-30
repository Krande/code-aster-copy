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
    subroutine amumpp(option, nbsol, kxmps, ldist, type, &
                      impr, ifmump, eli2lg, rsolu, csolu, &
                      vcine, prepos, lpreco, lmhpc)
        integer(kind=4) :: option
        integer(kind=8) :: nbsol
        integer(kind=8) :: kxmps
        aster_logical :: ldist
        character(len=1) :: type
        character(len=14) :: impr
        integer(kind=8) :: ifmump
        aster_logical :: eli2lg
        real(kind=8) :: rsolu(*)
        complex(kind=8) :: csolu(*)
        character(len=19) :: vcine
        aster_logical :: prepos
        aster_logical :: lpreco
        aster_logical :: lmhpc
    end subroutine amumpp
end interface
