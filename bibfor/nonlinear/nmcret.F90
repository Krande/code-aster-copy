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
subroutine nmcret(sderro, errorCodeType, errorCodeVale)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcrel.h"
#include "asterfort/NonLinear_type.h"
!
    character(len=24), intent(in) :: sderro
    character(len=3), intent(in) :: errorCodeType
    integer(kind=8), intent(in) :: errorCodeVale
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (SD ERREUR)
!
! GENERE UN EVENEMENT A PARTIR D'UN CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
! In  sderro           : name of datastructure for events in algorithm
! IN  errorCodeType : TYPE DU CODE RETOUR
!             'LDC' - INTEG. DU COMPORTEMENT
!                -1 : PAS D'INTEGRATION DU COMPORTEMENT
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE L'INTEGRATION DE LA LDC
!                 2 : ERREUR SUR LA NON VERIF. DE CRITERES PHYSIQUES
!                 3 : SIZZ PAS NUL POUR C_PLAN DEBORST
!             'PIL' - PILOTAGE
!                -1 : PAS DE CALCUL DU PILOTAGE
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : PAS DE SOLUTION
!                 2 : BORNE ATTEINTE -> FIN DU CALCUL
!             'FAC' - FACTORISATION
!                -1 : PAS DE FACTORISATION
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : MATRICE SINGULIERE
!                 2 : ERREUR LORS DE LA FACTORISATION
!                 3 : ON NE SAIT PAS SI SINGULIERE
!             'RES' - RESOLUTION
!                -1 : PAS DE RESOLUTION
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : NOMBRE MAXIMUM D'ITERATIONS ATTEINT
!             'CTC' - CONTACT DISCRET
!                -1 : PAS DE CALCUL DU CONTACT DISCRET
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : NOMBRE MAXI D'ITERATIONS
!                 2 : MATRICE SINGULIERE
! IN  errorCodeVale   : VALEUR DU CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iEvent
    character(len=24) :: eventECONJv, eventECOVJv, eventENOMJv
    integer(kind=8), pointer :: eventECOV(:) => null()
    character(len=16), pointer :: eventENOM(:) => null()
    character(len=8), pointer :: eventECON(:) => null()
    character(len=9) :: eventName
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    if (errorCodeVale .eq. -1) goto 999

! - Access to datastructure
    eventENOMJv = sderro(1:19)//'.ENOM'
    eventECOVJv = sderro(1:19)//'.ECOV'
    eventECONJv = sderro(1:19)//'.ECON'
    call jeveuo(eventENOMJv, 'L', vk16=eventENOM)
    call jeveuo(eventECOVJv, 'L', vi=eventECOV)
    call jeveuo(eventECONJv, 'L', vk8=eventECON)

! - Look for event
    do iEvent = 1, ZEVEN
        eventName = eventENOM(iEvent) (1:9)
        if (eventECON(iEvent) .eq. errorCodeType) then
            if (eventECOV(iEvent) .eq. errorCodeVale) then
                call nmcrel(sderro, eventName, .true._1)
            else
                call nmcrel(sderro, eventName, .false._1)
            end if
        end if
    end do
!
999 continue
!
    call jedema()
end subroutine
