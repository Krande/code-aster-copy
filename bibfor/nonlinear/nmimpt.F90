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
subroutine nmimpt(nume_inst, sddisc, ds_print)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/dinins.h"
#include "asterfort/nmimcr.h"
#include "asterfort/nmimen.h"
#include "asterfort/nmimr0.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "asterfort/jeexin.h"
!
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: sddisc
    type(NL_DS_Print), intent(inout) :: ds_print
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Print management
!
! Print head and convergence table
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_inst        : index of current time step
! In  sddisc           : name of datastructure for time discretization
! IO  ds_print         : datastructure for printing parameters
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: curr_inst
    integer(kind=8) :: lenivo, i_exist
    character(len=16) :: metlis
    character(len=24) :: sddisc_linf
!
! --------------------------------------------------------------------------------------------------
!
    sddisc_linf = sddisc(1:19)//'.LINF'
    call jeexin(sddisc_linf, i_exist)
!
! - Get level
!
    lenivo = 999
    if (i_exist .eq. 0) then
        lenivo = 0
    elseif (nume_inst .eq. 0) then
        lenivo = 0
    else
        call utdidt('L', sddisc, 'LIST', 'METHODE', valk_=metlis)
        if (metlis .eq. 'MANUEL') then
            lenivo = dinins(sddisc, nume_inst)
        else if (metlis .eq. 'AUTO') then
            lenivo = 0
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
! - Get current time and set in row
!
    curr_inst = diinst(sddisc, nume_inst)
    call nmimcr(ds_print, 'INCR_INST', curr_inst, .true._1)
!
! - Set values are not affected on rows for time loop
!
    call nmimr0(ds_print, 'INST')
!
! - Print current time
!
    if (lenivo .eq. 0) then
        call utmess('I', 'MECANONLINE6_6', sr=curr_inst)
    else
        call utmess('I', 'MECANONLINE6_1', si=lenivo, sr=curr_inst)
    end if
!
! - Print head of convergence table
!
    call nmimen(ds_print)
!
end subroutine
