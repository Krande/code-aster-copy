! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1403
!
module NonLinearDyna_type
! ==================================================================================================
    use Damping_type
! ==================================================================================================
    implicit none
!
#include "asterf_types.h"
!
! --------------------------------------------------------------------------------------------------
!
! Non-Linear operator for dynamic
!
! Define types
!
! --------------------------------------------------------------------------------------------------
!
! - Type: Damping parameters
    type NLDYNA_DAMPING
        aster_logical :: lDampRayleigh = ASTER_FALSE
        aster_logical :: lDampRayleighTang = ASTER_FALSE
        aster_logical :: lDampContact = ASTER_FALSE
        aster_logical :: lDampFEModel = ASTER_FALSE
        aster_logical :: lDampDiscret = ASTER_FALSE
        aster_logical :: lDampModal = ASTER_FALSE
        aster_logical :: hasMatrDamp = ASTER_FALSE
! - For modal damping
        type(MODAL_DAMPING) :: modalDamping
    end type NLDYNA_DAMPING
!
end module NonLinearDyna_type
