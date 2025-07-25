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

subroutine nmtime(ds_measure, operation_, device_type_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/nmrtim.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/GetDevice.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=*), intent(in) :: operation_
    character(len=*), intent(in) :: device_type_
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Measure and statistic management
!
! Timer management for device
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_measure       : datastructure for measure and statistics management
! In  operation        : operation to do Init, Launch or stop timer for current device
! In  device_type      : type of current device
!
! --------------------------------------------------------------------------------------------------
!
    type(NL_DS_Device) :: device
    character(len=9) :: timer_type
    character(len=10) :: device_type
    integer(kind=8) :: i_timer, timer_indx, nb_timer
    type(NL_DS_Timer) :: timer
    character(len=24) :: operation, cpu_name
    real(kind=8) :: time, list_time(7)
!
! --------------------------------------------------------------------------------------------------
!
    operation = operation_
    list_time(1:7) = 0.d0
!
! - Get current device
!
    call GetDevice(ds_measure, device_type_, device)
    timer_type = device%timer_name
    device_type = device%type
!
! - Find timer
!
    nb_timer = ds_measure%nb_timer
    timer_indx = 0
    do i_timer = 1, nb_timer
        if (ds_measure%timer(i_timer)%type .eq. timer_type) then
            ASSERT(timer_indx .eq. 0)
            timer_indx = i_timer
        end if
    end do
!
! - Get current timer
!
    ASSERT(timer_indx .ne. 0)
    timer = ds_measure%timer(timer_indx)
    cpu_name = timer%cpu_name
!
! - Operations
!
    if (operation .eq. 'Init') then
        call uttcpu(cpu_name, 'INIT', ' ')
        timer%time_init = 0.d0
    elseif (operation .eq. 'Launch') then
        call uttcpu(cpu_name, 'DEBUT', ' ')
    elseif (operation .eq. 'Stop') then
        call uttcpu(cpu_name, 'FIN', ' ')
        call uttcpr(cpu_name, 7, list_time)
        time = list_time(7)-timer%time_init
        timer%time_init = list_time(7)
        if (timer_type .eq. 'Store') then
            ds_measure%store_mean_time = list_time(4)
        end if
        if (timer_type .eq. 'Time_Step') then
            ds_measure%step_mean_time = list_time(4)
            ds_measure%step_remain_time = list_time(1)
        end if
        if (timer_type .eq. 'Newt_Iter') then
            ds_measure%iter_mean_time = list_time(4)
            ds_measure%iter_remain_time = list_time(1)
        end if
        call nmrtim(ds_measure, device_type, time)
    else
        ASSERT(.false.)
    end if
!
! - Save timer
!
    ds_measure%timer(timer_indx) = timer
!
end subroutine
