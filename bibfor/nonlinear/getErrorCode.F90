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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine getErrorCode(codret, ldccvg)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nmiret.h"
!
    character(len=19), intent(in) :: codret
    integer(kind=8), intent(out) :: ldccvg
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: tabret(0:10)
!
! --------------------------------------------------------------------------------------------------
!
    ldccvg = 0
    call nmiret(codret, tabret)
    if (tabret(0)) then
        if (tabret(4)) then
            ldccvg = 4
        else if (tabret(3)) then
            ldccvg = 3
        else if (tabret(2)) then
            ldccvg = 2
        else
            ldccvg = 1
        end if
        if (tabret(1)) then
            ldccvg = 1
        end if
    end if
    !
end subroutine
