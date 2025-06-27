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
subroutine setTimeListProgressBar(sddisc, nume_inst, final_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/getTimeListBounds.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: nume_inst
    aster_logical, optional, intent(in) :: final_
!
! --------------------------------------------------------------------------------------------------
!
! Display progress bar for time step list
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current inst step
! In  final            : flag for final time step
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: perc, numeStore
    real(kind=8) :: time_curr, t_ini, t_end, timeStore
    character(len=24) :: sdarchAinfJv
    integer(kind=8), pointer :: sdarchAinf(:) => null()
    character(len=19) :: sdarch
    character(len=24) :: sdarchListJv
    real(kind=8), pointer :: sdarchList(:) => null()
    integer(kind=8) :: iret, stateStoring
    aster_logical :: hasStoringList, lStoringInitState, lStoringSomething
!
! --------------------------------------------------------------------------------------------------
!
    time_curr = diinst(sddisc, nume_inst)
    sdarch = sddisc(1:14)//'.ARCH'
    sdarchAinfJv = sdarch(1:19)//'.AINF'
    sdarchListJv = sdarch(1:19)//'.LIST'
    call jeveuo(sdarchAinfJv, 'L', vi=sdarchAinf)

    call jeexin(sdarchListJv, iret)
    hasStoringList = iret .gt. 0

! - Compute percentage
    call getTimeListBounds(sddisc, t_ini, t_end)
    if (present(final_)) then
        ASSERT(final_)
        perc = 100.
    else
        perc = int(100.d0*(time_curr-t_ini)/(t_end-t_ini))
    end if

! - Get time
    lStoringInitState = ASTER_FALSE
    lStoringSomething = ASTER_TRUE
    if (hasStoringList) then
        call jeveuo(sdarchListJv, 'L', vr=sdarchList)
        stateStoring = sdarchAinf(4)
        lStoringInitState = ASTER_FALSE
        lStoringSomething = ASTER_TRUE
        if (stateStoring .eq. -1) then
            lStoringSomething = ASTER_FALSE
            timeStore = 0.d0
            numeStore = -1

        elseif (stateStoring .eq. 0) then
            lStoringInitState = ASTER_TRUE
            timeStore = 0.d0
            numeStore = 0

        else
            lStoringSomething = ASTER_TRUE
            numeStore = stateStoring
            timeStore = sdarchList(numeStore)

        end if
    else
        lStoringInitState = ASTER_FALSE
        lStoringSomething = ASTER_TRUE
        timeStore = time_curr
        numeStore = nume_inst
    end if

! - Cas où l'on force l'archivage du dernier instant: il n'est pas dans la liste !
    if (present(final_)) then
        lStoringInitState = ASTER_FALSE
        lStoringSomething = ASTER_TRUE
        timeStore = time_curr
    end if

    if (lStoringInitState) then
        call utmess('I', 'PROGRESS_4', si=perc, sr=time_curr)

    elseif (lStoringSomething) then
        call utmess('I', 'PROGRESS_1', ni=2, vali=[perc, numeStore], &
                    nr=2, valr=[time_curr, timeStore])

    else
        call utmess('I', 'PROGRESS_3', si=perc, sr=time_curr)

    end if
!
end subroutine
