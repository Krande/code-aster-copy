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
    subroutine aceinc(noma, nomo, ntyele, nbocc, ivr, &
                      locaco, locagb, locamb, zjdlm, lmax, ier)
        character(len=8) :: noma
        character(len=8) :: nomo
        integer(kind=8) :: ntyele(*)
        integer(kind=8) :: nbocc(*)
        integer(kind=8) :: ivr(*)
        aster_logical :: locaco
        aster_logical :: locagb
        aster_logical :: locamb
        integer(kind=8) :: zjdlm(*)
        integer(kind=8) :: lmax
        integer(kind=8) :: ier
    end subroutine aceinc
end interface
