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

subroutine nmtima(ds_measure, timer_type_, vali)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Measure), intent(in) :: ds_measure
    character(len=*), intent(in) :: timer_type_
    integer(kind=8), intent(out) :: vali
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Measure and statistic management
!
! Evaluate remaining type
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_measure       : datastructure for measure and statistics management
! In  timer_type       : type of current timer
! Out vali             : 0 - enough time
!                        1 - not enough time
!
! --------------------------------------------------------------------------------------------------
!
    character(len=9) :: timer_type
    integer(kind=8) :: i_timer, timer_indx, nb_timer
    type(NL_DS_Timer) :: timer
    real(kind=8) :: remaining_time, store_mean_time, iter_mean_time, step_mean_time
!
! --------------------------------------------------------------------------------------------------
!
    nb_timer = ds_measure%nb_timer
    timer_indx = 0
    vali = 0
    timer_type = timer_type_
!
! - Find timer
!
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
!
! - Get mean times
!
    store_mean_time = ds_measure%store_mean_time
    iter_mean_time = ds_measure%iter_mean_time
    step_mean_time = ds_measure%step_mean_time
!
! - Enough time ?
!
    if (timer_type .eq. 'Newt_Iter') then
        remaining_time = ds_measure%iter_remain_time
        if ((2.d0*iter_mean_time) .le. (0.95d0*remaining_time-store_mean_time)) then
            vali = 0
        else
            vali = 1
        end if
    else if (timer_type .eq. 'Time_Step') then
        remaining_time = ds_measure%step_remain_time
        if (step_mean_time .le. 0.90d0*remaining_time) then
            vali = 0
        else
            vali = 1
        end if
    else
        ASSERT(.false.)
    end if
!
end subroutine
