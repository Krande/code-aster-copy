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
subroutine nmltev(sderro, typevt, loopName, eventFlag)
!
    implicit none
    !
#include "asterf_types.h"
#include "asterfort/asmpi_any.h"
#include "asterfort/jeveuo.h"
#include "asterfort/NonLinear_type.h"
!
    character(len=24) :: sderro
    character(len=4) :: typevt
    character(len=4) :: loopName
    aster_logical :: eventFlag
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (SD ERREUR)
!
! DIT SI UN EVENEMENT DE TYPE DONNE EST ACTIVE
!
! --------------------------------------------------------------------------------------------------
!
! IN  SDERRO : SD GESTION DES ERREURS
! IN  TYPEVT : TYPE EVENEMENT
!               'EVEN' - EVENEMENT SIMPLE
!               'ERRI' - EVENEMENT DE TYPE ERREUR IMMEDIATE
!               'ERRC' - EVENEMENT DE TYPE ERREUR A CONVERGENCE
!               'CONV' - EVENEMENT POUR DETERMINER LA CONVERGENCE
! IN  NOMBCL : NOM DE LA BOUCLE
!               'RESI' - RESIDUS D'EQUILIBRE
!               'NEWT' - BOUCLE DE NEWTON
!               'FIXE' - BOUCLE DE POINT FIXE
!               'INST' - BOUCLE SUR LES PAS DE TEMPS
! OUT LEVENT : .TRUE. SI AU MOINS UN EVENT EST ACTIVE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iEvent, eventState
    character(len=9) :: eventLevel
    character(len=24) :: eventEACTJv, eventENIVJv
    integer(kind=8), pointer :: eventEACT(:) => null()
    character(len=16), pointer :: eventENIV(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    eventFlag = .false.

! - Access to datastructure
    eventEACTJv = sderro(1:19)//'.EACT'
    eventENIVJv = sderro(1:19)//'.ENIV'
    call jeveuo(eventEACTJv, 'L', vi=eventEACT)
    call jeveuo(eventENIVJv, 'L', vk16=eventENIV)
!
! --- AU MOINS UN EVENEMENT DE CE NIVEAU D'ERREUR EST ACTIVE ?
!
    do iEvent = 1, ZEVEN
        eventState = eventEACT(iEvent)
        eventLevel = eventENIV(iEvent) (1:9)
        if (eventLevel(1:4) .eq. typevt) then
            if (typevt .eq. 'EVEN') then
                if (eventState .eq. EVENT_IS_ACTIVE) then
                    eventFlag = .true.
                end if
            else
                if (eventLevel(6:9) .eq. loopName) then
                    if (eventState .eq. EVENT_IS_ACTIVE) then
                        eventFlag = .true.
                    end if
                end if
            end if
        end if
    end do

! - Share error for HPC
    eventFlag = asmpi_any(eventFlag, ASTER_TRUE)
!
end subroutine
