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
#include "asterf_types.h"
!
interface
    subroutine tanbul(option, ndim, g, mate, compor,&
                      resi, mini, alpha, dsbdep, trepst)
        integer(kind=8) :: ndim
        character(len=16) :: option
        integer(kind=8) :: g
        integer(kind=8) :: mate
        character(len=16) :: compor
        aster_logical :: resi
        aster_logical :: mini
        real(kind=8) :: alpha
        real(kind=8) :: dsbdep(2*ndim, 2*ndim)
        real(kind=8) :: trepst
    end subroutine tanbul
end interface
