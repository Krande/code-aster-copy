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

subroutine nmtimr(ds_measure, device_type_, phasis, time)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/GetDevice.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Measure), intent(in) :: ds_measure
    character(len=*), intent(in) :: device_type_
    character(len=1), intent(in) :: phasis
    real(kind=8), intent(out) :: time
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Measure and statistic management
!
! Get time
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_measure       : datastructure for measure and statistics management
! In  device_type      : type of current device
! In  phasis           : phasis (time step, Newton iteration, all computation)
! Out time             : time
!
! --------------------------------------------------------------------------------------------------
!
    type(NL_DS_Device) :: device
    real(kind=8) :: time_iter, time_step, time_comp
!
! --------------------------------------------------------------------------------------------------
!
!
! - Get current device
!
    call GetDevice(ds_measure, device_type_, device)
!
! - Get times
!
    time_iter = device%time_iter
    time_step = device%time_step
    time_comp = device%time_comp
    if (phasis .eq. 'T') then
        time = time_comp
    elseif (phasis .eq. 'P') then
        time = time_step
    elseif (phasis .eq. 'N') then
        time = time_iter
    else
        ASSERT(.false.)
    end if
!
end subroutine
