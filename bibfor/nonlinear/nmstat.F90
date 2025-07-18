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
subroutine nmstat(phasis, ds_measure, ds_print, sddisc, nume_inst, sderro)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/diinst.h"
#include "asterfort/impmem.h"
#include "asterfort/GetDevice.h"
#include "asterfort/nonlinDSPrintTableLine.h"
#include "asterfort/nmrini.h"
#include "asterfort/nmrtim.h"
#include "asterfort/nmtimr.h"
#include "asterfort/nmstat_mess.h"
#include "asterfort/nmstat_table.h"
#include "asterfort/nmstat_vale.h"
!
    character(len=1), intent(in) :: phasis
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_Print), intent(in) :: ds_print
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: nume_inst
    character(len=24), intent(in) :: sderro
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Measure and statistics management
!
! Update statistics
!
! --------------------------------------------------------------------------------------------------
!
! In  phasis           : phasis
!                          'P' current step time
!                          'T' on all transient
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_print         : datastructure for printing parameters
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current step time
! In  sderro           : datastructure for errors during algorithm
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_print
    real(kind=8) :: time_other, time_curr
    real(kind=8) :: time_time_step
    real(kind=8) :: time_sub, time
    integer(kind=8) :: i_device, nb_device
    type(NL_DS_Device) :: device
    character(len=10) :: device_type
!
! --------------------------------------------------------------------------------------------------
!
    time_other = 0.d0
    time_time_step = 0.d0
    time_sub = 0.d0
    time = 0.d0
!
! - Current time
!
    time_curr = diinst(sddisc, nume_inst)
!
! - Print for this step time ?
!
    l_print = ds_print%l_print
!
! - Time for other operations (not measured)
!
    nb_device = ds_measure%nb_device
    if (phasis .eq. 'P') then
        call nmtimr(ds_measure, 'Time_Step', phasis, time_time_step)
        time_sub = 0.d0
        do i_device = 1, nb_device
            device = ds_measure%device(i_device)
            device_type = device%type
            call GetDevice(ds_measure, device_type, device)
            if (device%time_indi_step .ne. 0 .and. device_type .ne. 'Time_Step') then
                call nmtimr(ds_measure, device_type, phasis, time)
                time_sub = time_sub+time
            end if
        end do
        time_other = time_time_step-time_sub
        if (time_other .le. 0.d0) then
            time_other = 0.d0
        end if
        call nmrtim(ds_measure, 'Other', time_other)
    end if
!
! - Save values in columns
!
    if (phasis .eq. 'P') then
        call nmstat_vale(ds_measure, time_curr, sderro)
    end if
!
! - Print at end of current step time
!
    if ((phasis .eq. 'P') .and. l_print) then
        call nmstat_mess(ds_measure, phasis)
        call impmem()
    end if
!
! - Save in table
!
    if ((phasis .eq. 'P') .and. ds_measure%l_table) then
        call nmstat_table(ds_measure)
    end if
!
! - Save in file
!
    if ((phasis .eq. 'P') .and. ds_measure%table%l_csv) then
        call nonlinDSPrintTableLine(ds_measure%table, ',', ds_measure%table%unit_csv)
    end if
!
! - Print at end of computation
!
    if (phasis .eq. 'T') then
        call nmstat_mess(ds_measure, phasis)
    end if
!
! - Reset times and counters
!
    call nmrini(ds_measure, phasis)
!
end subroutine
