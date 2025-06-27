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
subroutine nmeraz(sderro, eventType)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/NonLinear_type.h"
#include "jeveux.h"
!
    character(len=24), intent(in) :: sderro
    character(len=4), intent(in) :: eventType
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (SD ERREUR)
!
! REMISE A ZERO DES EVENEMENTS
!
! --------------------------------------------------------------------------------------------------
!
! IN  SDERRO : SD GESTION DES ERREURS
! IN  TYPEVT : TYPE DE L'EVENEMENT
!              'TOUS' - TOUS LES EVENEMENTS
!              'EVEN' - EVENEMENT SIMPLE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iEvent, iret
    character(len=16) :: eventLevel
    character(len=24) :: eventEACTJv, eventENIVJv
    integer(kind=8), pointer :: eventEACT(:) => null()
    character(len=16), pointer :: eventENIV(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    eventEACTJv = sderro(1:19)//'.ENOM'
    call jeexin(eventEACTJv, iret)
    if (iret .eq. 0) goto 99

! - Access to datastructure
    eventEACTJv = sderro(1:19)//'.EACT'
    eventENIVJv = sderro(1:19)//'.ENIV'
    call jeveuo(eventEACTJv, 'L', vi=eventEACT)
    call jeveuo(eventENIVJv, 'L', vk16=eventENIV)

! - EVENEMENTS DESACTIVES
    do iEvent = 1, ZEVEN
        eventLevel = eventENIV(iEvent) (1:9)
        if (eventType .eq. 'TOUS') then
            eventEACT(iEvent) = EVENT_IS_INACTIVE
        else if (eventType .eq. 'EVEN') then
            if (eventLevel .eq. 'EVEN') then
                eventEACT(iEvent) = EVENT_IS_INACTIVE
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
99  continue
!
    call jedema()
end subroutine
