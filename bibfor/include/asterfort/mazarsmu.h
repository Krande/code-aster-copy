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
    subroutine mazarsmu(option, epsela, deps, dimloc, mazars, varm, varp, sigp, dsidep)
        character(len=16)   :: option
        real(kind=8)        :: epsela(6)
        real(kind=8)        :: deps(6)
        integer(kind=8)             :: dimloc
        real(kind=8)        :: mazars(*)
        real(kind=8)        :: varm(*)
        real(kind=8)        :: varp(*)
        real(kind=8)        :: sigp(*)
        real(kind=8)        :: dsidep(6, 6)
    end subroutine mazarsmu
end interface
