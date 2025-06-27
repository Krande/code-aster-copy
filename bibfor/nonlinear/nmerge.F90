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
subroutine nmerge(sderro, eventToTest, eventFlag)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/NonLinear_type.h"
!
    character(len=24), intent(in) :: sderro
    character(len=9), intent(in) :: eventToTest
    aster_logical, intent(out) :: eventFlag
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (SD ERREUR)
!
! DIT SI UN EVENEMENT EST DECLENCHE
!
! --------------------------------------------------------------------------------------------------
!
! In  sderro           : name of datastructure for events in algorithm
! IN  NOMEVT : NOM DE L'EVENEMENT (VOIR LA LISTE DANS NMCRER)
! OUT LACTIV : .TRUE. SI EVENEMENT ACTIVE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iEvent, eventState
    character(len=24) :: eventENOMJv, eventEACTJv
    integer(kind=8), pointer :: eventEACT(:) => null()
    character(len=16), pointer :: eventENOM(:) => null()
    character(len=16) :: eventName
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to datastructure
    eventENOMJv = sderro(1:19)//'.ENOM'
    eventEACTJv = sderro(1:19)//'.EACT'
    call jeveuo(eventENOMJv, 'L', vk16=eventENOM)
    call jeveuo(eventEACTJv, 'L', vi=eventEACT)

! - Test
    eventFlag = ASTER_FALSE
    do iEvent = 1, ZEVEN
        eventName = eventENOM(iEvent)
        if (eventName .eq. eventToTest) then
            eventState = eventEACT(iEvent)
            if (eventState .eq. EVENT_IS_ACTIVE) then
                eventFlag = ASTER_TRUE
            end if
        end if
    end do
!
    call jedema()
end subroutine
