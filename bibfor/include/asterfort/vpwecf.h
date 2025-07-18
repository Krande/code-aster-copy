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
    subroutine vpwecf(option, typres, nfreq, mxfreq, resufi,&
                      resufr, resufk, lamor, ktyp, lns)
        integer(kind=8) :: mxfreq
        character(len=*) :: option
        character(len=*) :: typres
        integer(kind=8) :: nfreq
        integer(kind=8) :: resufi(mxfreq, *)
        real(kind=8) :: resufr(mxfreq, *)
        character(len=*) :: resufk(mxfreq, *)
        integer(kind=8) :: lamor
        character(len=1) :: ktyp
        aster_logical :: lns
    end subroutine vpwecf
end interface
