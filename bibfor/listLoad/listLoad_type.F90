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
! ==================================================================================================
!
! Types for management of list of loads
!
! ==================================================================================================
!
module listLoad_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
! ==================================================================================================
! Type: parameters to construct list of loads
! ==================================================================================================
    type ListLoad_Prep
! ----- Model from command
        character(len=8) :: model = " "
! ----- Command STAT_NON_LINE
        aster_logical :: staticOperator = ASTER_FALSE
! ----- Type (real or complex) for function multiplier
        aster_logical :: funcIsCplx = ASTER_FALSE
! ----- Continuation method (PILOTAGE)
        aster_logical :: lHasPilo = ASTER_FALSE
    end type ListLoad_Prep
!===================================================================================================
    public :: ListLoad_Prep
contains
!===================================================================================================
end module listLoad_type
