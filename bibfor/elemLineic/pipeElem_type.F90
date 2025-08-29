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
! Types for pipe elements
!
! ==================================================================================================
!
module pipeElem_type
! ==================================================================================================
    use beamElem_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/pipeElem_type.h"
! ==================================================================================================
! Type: properties of pipe element
! ==================================================================================================
    type pipeElem_Prop
        integer(kind=8) :: nbNode = 0
! ----- Fourier modes
        integer(kind=8) :: nbFourier = 0
! ----- Properties of beam
        type(beamElem_Prop) :: beamElem
! ----- Global properties
        aster_logical :: lModiMetric = ASTER_FALSE
        integer(kind=8) :: pipeType = PIPE_TYPE_UNDEF
    end type pipeElem_Prop
!===================================================================================================
    public :: pipeElem_Prop
contains
!===================================================================================================
end module pipeElem_type
