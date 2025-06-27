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
subroutine nmecev(sderro, acces, failType, actionType)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/NonLinear_type.h"
!
    character(len=24), intent(in) :: sderro
    character(len=1), intent(in) :: acces
    integer(kind=8), intent(inout) :: failType
    integer(kind=8), intent(inout) :: actionType
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! ECHEC DU TRAITEMENT D'UNE ACTION - SAUVEGARDE/LECTURE POUR INFO
!
! --------------------------------------------------------------------------------------------------
!
! IN  SDERRO : SD ERREUR
! IN  ACCES  : TYPE ACCES 'E' OU 'L'
! I/O NOMEVD : NOM DE L'EVENEMENT
! I/O ACTION : NOM DE L'ACTION
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: eventEEVTJv
    integer(kind=8), pointer :: eventEEVT(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    eventEEVTJv = sderro(1:19)//'.EEVT'
    call jeveuo(eventEEVTJv, 'E', vi=eventEEVT)
    if (acces .eq. 'E') then
        eventEEVT(1) = failType
        eventEEVT(2) = actionType
    else if (acces .eq. 'L') then
        failType = eventEEVT(1)
        actionType = eventEEVT(2)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
