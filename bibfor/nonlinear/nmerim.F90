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
subroutine nmerim(sderro)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=24), intent(in) :: sderro
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (SD ERREUR)
!
! EMISSION MESSAGE ERRREUR
!
! --------------------------------------------------------------------------------------------------
!
! IN  SDERRO : SD ERREUR
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iEvent
    integer(kind=8) :: eventState
    character(len=9) :: eventLevel
    character(len=24) :: eventMesg
    character(len=24) :: eventEACTJv, eventENIVJv, eventEMSGJv
    integer(kind=8), pointer :: eventEACT(:) => null()
    character(len=16), pointer :: eventENIV(:) => null()
    character(len=24), pointer :: eventEMSG(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Access to datastructure
    eventEACTJv = sderro(1:19)//'.EACT'
    eventENIVJv = sderro(1:19)//'.ENIV'
    eventEMSGJv = sderro(1:19)//'.EMSG'
    call jeveuo(eventEACTJv, 'L', vi=eventEACT)
    call jeveuo(eventENIVJv, 'L', vk16=eventENIV)
    call jeveuo(eventEMSGJv, 'L', vk24=eventEMSG)

! - Print error messages
    do iEvent = 1, ZEVEN
        eventState = eventEACT(iEvent)
        eventLevel = eventENIV(iEvent) (1:9)
        eventMesg = eventEMSG(iEvent)
        if ((eventLevel(1:3) .eq. 'ERR') .and. (eventState .eq. EVENT_IS_ACTIVE)) then
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
            else if (eventMesg .eq. 'MECANONLINE10_13') then
                call utmess('I', 'MECANONLINE10_13')
            else if (eventMesg .eq. 'MECANONLINE10_20') then
                call utmess('I', 'MECANONLINE10_20')
            else if (eventMesg .eq. 'MECANONLINE10_24') then
                call utmess('I', 'MECANONLINE10_24')
            else if (eventMesg .eq. 'MECANONLINE10_36') then
                call utmess('I', 'MECANONLINE10_36')
            else if (eventMesg .eq. 'MECANONLINE10_14') then
                call utmess('I', 'MECANONLINE10_14')
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
    end do
!
end subroutine
