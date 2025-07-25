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
subroutine nmstat_vale(ds_measure, time_curr, sderro)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmtimr.h"
#include "asterfort/nmrvai.h"
#include "asterfort/utgtme.h"
#include "asterfort/SetTableColumn.h"
#include "asterfort/NonLinear_type.h"
!
    type(NL_DS_Measure), intent(inout) :: ds_measure
    real(kind=8), intent(in) :: time_curr
    character(len=24), intent(in) :: sderro
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Measure and statistics management
!
! Update statistics in columns
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_measure       : datastructure for measure and statistics management
! In  time_curr        : current time
! In  sderro           : datastructure for errors during algorithm
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: vmpeak(1)
    integer(kind=8) :: iret
    integer(kind=8) :: nb_cols, nb_device
    integer(kind=8) :: i_col, i_device
    type(NL_DS_Table) :: table
    type(NL_DS_Column) :: column
    type(NL_DS_Device) :: device
    aster_logical :: l_vale_inte, l_vale_real
    integer(kind=8) :: count, iEvent, eventState
    character(len=10) :: device_type
    real(kind=8) :: time
    character(len=16) :: col_name, eventName, state
    character(len=24) :: eventENOMJv, eventEACTJv
    integer(kind=8), pointer :: eventEACT(:) => null()
    character(len=16), pointer :: eventENOM(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    table = ds_measure%table
    nb_cols = table%nb_cols
    nb_device = ds_measure%nb_device

! - Access to datastructure
    eventENOMJv = sderro(1:19)//'.ENOM'
    eventEACTJv = sderro(1:19)//'.EACT'
    call jeveuo(eventENOMJv, 'L', vk16=eventENOM)
    call jeveuo(eventEACTJv, 'L', vi=eventEACT)

! - State of step
    state = 'CONV'

    do iEvent = 1, ZEVEN
        eventName = eventENOM(iEvent)
        eventState = eventEACT(iEvent)
        if (eventState .eq. EVENT_IS_ACTIVE) then
            state = eventName
            exit
        end if
    end do

! - Get memory
    call utgtme(1, 'VMPEAK  ', vmpeak, iret)

! - Set list of values in columns
    do i_col = 1, nb_cols
        column = table%cols(i_col)
        i_device = table%indx_vale(i_col)
        col_name = column%name
        if (i_device .eq. 0) then
            if (col_name .eq. 'INST') then
                call SetTableColumn(table, 'INST', flag_affe_=.true._1, valer_=time_curr)
            elseif (col_name .eq. 'State') then
                call SetTableColumn(table, 'State', flag_affe_=.true._1, valek_=state)
            elseif (col_name .eq. 'Memory') then
                call SetTableColumn(table, 'Memory', flag_affe_=.true._1, &
                                    valei_=nint(vmpeak(1)))
            end if
        else
            device = ds_measure%device(i_device)
            device_type = device%type
            l_vale_inte = column%l_vale_inte
            l_vale_real = column%l_vale_real
            if (l_vale_real) then
                call nmtimr(ds_measure, device_type, 'P', time)
                ASSERT(col_name(1:5) .eq. 'Time_')
                call SetTableColumn(table, col_name, flag_affe_=.true._1, valer_=time)
            end if
            if (l_vale_inte) then
                call nmrvai(ds_measure, device_type, 'P', output_count=count)
                ASSERT(col_name(1:6) .eq. 'Count_')
                call SetTableColumn(table, col_name, flag_affe_=.true._1, valei_=count)
            end if
        end if
    end do
!
    ds_measure%table = table
!
end subroutine
