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
! Types for beam elements
!
! ==================================================================================================
!
module beamElem_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    private
#include "asterf_types.h"
! ==================================================================================================
! Type: properties of pipe section
! ==================================================================================================
    type sectPipe_Prop
        real(kind=8) :: thickness = 0.d0
        real(kind=8) :: radiusExt = 0.d0
        real(kind=8) :: radiusInt = 0.d0
        real(kind=8) :: radiusMoy = 0.d0
        real(kind=8) :: area = 0.d0
    end type sectPipe_Prop
!===================================================================================================
    public :: sectPipe_Prop
contains
!===================================================================================================
end module beamElem_type
