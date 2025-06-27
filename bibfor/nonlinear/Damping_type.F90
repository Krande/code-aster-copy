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
module Damping_type
!
    implicit none
!
#include "asterf_types.h"
!
! --------------------------------------------------------------------------------------------------
!
! Damping
!
! Define types
!
! --------------------------------------------------------------------------------------------------
!

! - Type: modal damping parameters
    type MODAL_DAMPING
! - Flag for update speed
        aster_logical :: lReacVite = ASTER_FALSE
! - Modes for damping
        character(len=8) :: dampMode = " "
        integer(kind=8) :: nbMode = 0
! - Values of damping
        integer(kind=8) :: nbDampVale = 0
        character(len=24) :: jvListDamp = " "
! - Name of datastructure to save parameters
        character(len=24) :: jvDataDamp = " "
! - Flag for Debug
        aster_logical :: debug = ASTER_TRUE
    end type MODAL_DAMPING
!
end module Damping_type
