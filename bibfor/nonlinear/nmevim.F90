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
subroutine nmevim(ds_print, sddisc, sderro, loop_name)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmacto.h"
#include "asterfort/nmerge.h"
#include "asterfort/nmimpx.h"
#include "asterfort/nmlecv.h"
#include "asterfort/nmltev.h"
#include "asterfort/utmess.h"
#include "asterfort/getFailEvent.h"
#include "asterfort/NonLinear_type.h"
!
    type(NL_DS_Print), intent(in) :: ds_print
    character(len=24), intent(in) :: sderro
    character(len=19), intent(in) :: sddisc
    character(len=4), intent(in) :: loop_name
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Event management
!
! Print event messages
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_print         : datastructure for printing parameters
! In  sddisc           : datastructure for time discretization TEMPORELLE
! In  sderro           : name of datastructure for error management (events)
! In  loop_name        : name of loop
!                         'NEWT' - Newton loop
!                         'FIXE' - Fixed points loop
!                         'INST' - Step time loop
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lacti, cvbouc, lerrei, l_sep_line, lldcbo
    integer(kind=8) :: iFailActi
    integer(kind=8) :: iEvent
    character(len=24) :: eventEACTJv, eventENIVJv, eventEMSGJv
    integer(kind=8), pointer :: eventEACT(:) => null()
    character(len=16), pointer :: eventENIV(:) => null()
    character(len=24), pointer :: eventEMSG(:) => null()
    integer(kind=8) :: eventState, eventType
    character(len=9) :: eventLevel
    character(len=24) :: eventMesg
!
! --------------------------------------------------------------------------------------------------
!
    call nmlecv(sderro, loop_name, cvbouc)
    call nmltev(sderro, 'ERRI', loop_name, lerrei)
    call nmerge(sderro, 'INTE_BORN', lldcbo)

! - Separator line to print ?
    l_sep_line = (.not. cvbouc .and. .not. lerrei .and. .not. lldcbo)

! - Access to datastructure
    eventEACTJv = sderro(1:19)//'.EACT'
    eventENIVJv = sderro(1:19)//'.ENIV'
    eventEMSGJv = sderro(1:19)//'.EMSG'
    call jeveuo(eventEACTJv, 'L', vi=eventEACT)
    call jeveuo(eventENIVJv, 'L', vk16=eventENIV)
    call jeveuo(eventEMSGJv, 'L', vk24=eventEMSG)

! - Print event messages - Algorithm
    do iEvent = 1, ZEVEN
        eventState = eventEACT(iEvent)
        eventLevel = eventENIV(iEvent) (1:9)
        eventMesg = eventEMSG(iEvent)
        if (eventLevel(1:4) .eq. 'EVEN' .and. (eventState .eq. EVENT_IS_ACTIVE)) then
            if (eventMesg .ne. ' ') then
                if (l_sep_line) then
                    call nmimpx(ds_print)
                end if
                if (eventMesg .eq. 'MECANONLINE10_1') then
                    call utmess('I', 'MECANONLINE10_1')
                else if (eventMesg .eq. 'MECANONLINE10_2') then
                    call utmess('I', 'MECANONLINE10_2')
                else if (eventMesg .eq. 'MECANONLINE10_3') then
                    call utmess('I', 'MECANONLINE10_3')
                else if (eventMesg .eq. 'MECANONLINE10_4') then
                    call utmess('I', 'MECANONLINE10_4')
                else if (eventMesg .eq. 'MECANONLINE10_5') then
                    call utmess('I', 'MECANONLINE10_5')
                else if (eventMesg .eq. 'MECANONLINE10_6') then
                    call utmess('I', 'MECANONLINE10_6')
                else if (eventMesg .eq. 'MECANONLINE10_7') then
                    call utmess('I', 'MECANONLINE10_7')
                else if (eventMesg .eq. 'MECANONLINE10_8') then
                    call utmess('I', 'MECANONLINE10_8')
                else if (eventMesg .eq. 'MECANONLINE10_9') then
                    call utmess('I', 'MECANONLINE10_9')
                else if (eventMesg .eq. 'MECANONLINE10_10') then
                    call utmess('I', 'MECANONLINE10_10')
                else if (eventMesg .eq. 'MECANONLINE10_11') then
                    call utmess('I', 'MECANONLINE10_11')
                else if (eventMesg .eq. 'MECANONLINE10_12') then
                    call utmess('I', 'MECANONLINE10_12')
                else if (eventMesg .eq. 'MECANONLINE10_14') then
                    call utmess('I', 'MECANONLINE10_14')
                else if (eventMesg .eq. 'MECANONLINE10_20') then
                    call utmess('I', 'MECANONLINE10_20')
                else if (eventMesg .eq. 'MECANONLINE10_24') then
                    call utmess('I', 'MECANONLINE10_24')
                else if (eventMesg .eq. 'MECANONLINE10_26') then
                    call utmess('I', 'MECANONLINE10_26')
                else if (eventMesg .eq. 'MECANONLINE10_25') then
                    if (cvbouc .and. loop_name .eq. 'NEWT') then
                        call utmess('A', 'MECANONLINE10_25')
                    end if
                else
                    ASSERT(.false.)
                end if
            end if
        end if
    end do

! - Print event messages - User
    call nmacto(sddisc, iFailActi)
    lacti = iFailActi .gt. 0
    if (lacti) then
! ----- Get event type
        call getFailEvent(sddisc, iFailActi, eventType)
        if (eventType .eq. FAIL_EVT_INTERPENE) then
            if (l_sep_line) then
                call nmimpx(ds_print)
            end if
            call utmess('I', 'MECANONLINE10_22')
        else if (eventType .eq. FAIL_EVT_DIVE_RESI) then
            if (l_sep_line) then
                call nmimpx(ds_print)
            end if
            call utmess('I', 'MECANONLINE10_23')
        else if (eventType .eq. FAIL_EVT_RESI_MAXI) then
            if (l_sep_line) then
                call nmimpx(ds_print)
            end if
            call utmess('I', 'MECANONLINE10_26')
        else if (eventType .eq. FAIL_EVT_INCR_QUANT) then
            if (l_sep_line) then
                call nmimpx(ds_print)
            end if
            call utmess('I', 'MECANONLINE10_24')
        end if
    end if
!
end subroutine
