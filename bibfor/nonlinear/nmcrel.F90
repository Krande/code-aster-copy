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
subroutine nmcrel(sderro, eventName, eventFlag)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmeceb.h"
#include "asterfort/NonLinear_type.h"
!
    character(len=24), intent(in) :: sderro
    character(len=9), intent(in) :: eventName
    aster_logical, intent(in) :: eventFlag
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (SD ERREUR)
!
! ENREGISTREMENT D'UN EVENEMENT INTRINSEQUE
!
! --------------------------------------------------------------------------------------------------
!
! POUR LES EVENEMENTS A CODE RETOUR, IL FAUT D'ABORD TRANSFORMER LE
! CODE RETOUR EN EVENEMENT - UTILISER NMCRET
!
! In  sderro           : name of datastructure for events in algorithm
! In  eventName        : name of event
! In  eventFlag
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iEvent
    character(len=24) :: eventENOMJv, eventENIVJv, eventEACTJv
    integer(kind=8), pointer :: eventEACT(:) => null()
    character(len=16), pointer :: eventENOM(:) => null(), eventENIV(:) => null()
    character(len=16) :: eventLevel
    character(len=4) :: loopName
    integer(kind=8) :: eventIndx
!
! --------------------------------------------------------------------------------------------------
!
    eventIndx = 0

! - Access to datastructure
    eventENOMJv = sderro(1:19)//'.ENOM'
    eventENIVJv = sderro(1:19)//'.ENIV'
    eventEACTJv = sderro(1:19)//'.EACT'
    call jeveuo(eventENOMJv, 'L', vk16=eventENOM)
    call jeveuo(eventENIVJv, 'L', vk16=eventENIV)
    call jeveuo(eventEACTJv, 'E', vi=eventEACT)

! - look for event
    do iEvent = 1, ZEVEN
        if (eventENOM(iEvent) .eq. eventName) then
            eventIndx = iEvent
            goto 66
        end if
    end do
66  continue
    ASSERT(eventIndx .ne. 0)

! - (DES-)ACTIVATION DE L'EVENEMENT
    if (eventFlag) then
        eventEACT(eventIndx) = EVENT_IS_ACTIVE
    else
        eventEACT(eventIndx) = EVENT_IS_INACTIVE
    end if

! --- EVENEMENT DE TYPE ERREUR ACTIVE: ON CHANGE LE STATUT DE LA BOUCLE
    if (eventFlag) then
        eventLevel = eventENIV(eventIndx)
        if (eventLevel(1:5) .eq. 'ERRI_') then
            loopName = eventLevel(6:9)
            call nmeceb(sderro, loopName, 'ERRE')
            if (loopName .eq. 'CALC') then
                call nmeceb(sderro, 'NEWT', 'STOP')
                call nmeceb(sderro, 'FIXE', 'STOP')
                call nmeceb(sderro, 'INST', 'STOP')
            end if
        end if
    end if
!
end subroutine
