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
    subroutine rbph01(trange, nbcham, typea, itresu, nfonct,&
                      basemo, typref, typbas, tousno, multap, i_cham)
        character(len=19) :: trange
        integer(kind=8) :: nbcham
        character(len=16) :: typea(*)
        integer(kind=8) :: itresu(*)
        integer(kind=8) :: nfonct
        character(len=8) :: basemo
        character(len=19) :: typref(*)
        character(len=16) :: typbas(*)
        aster_logical :: tousno
        aster_logical :: multap
        integer(kind=8), dimension(8), intent(out), optional :: i_cham
    end subroutine rbph01
end interface
