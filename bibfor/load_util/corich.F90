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
subroutine corich(actionZ, fieldZ, ichin_, ichout_)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/wkvect.h"
!
    character(len=*), intent(in) :: actionZ, fieldZ
    integer(kind=8), optional, intent(in) :: ichin_
    integer(kind=8), optional, intent(out) :: ichout_
!
! --------------------------------------------------------------------------------------------------
!
!  BUT :   GERER UN EVENTUEL LIEN ENTRE UN CHAMP (RESUELEM OU CHAM_NO)
!          ET LE NUMERO DE LA CHARGE AUQUEL IL EST ASSOCIE POUR POUVOIR
!          LUI APPLIQUER LE BON FONC_MULT.
!   CE LIEN EST UTILISE PAR LA ROUTINE ASASVE (RESUELEM)
!   CE LIEN EST UTILISE PAR LA ROUTINE ASCOVA (CHAM_NO)
!
! IN  ACTION  K1 : / 'E' : ON ECRIT UN LIEN
!                  / 'L' : ON LIT UN LIEN ECRIT AU PREALABLE
!                  / 'S' : ON SUPPRIME UN LIEN ECRIT AU PREALABLE
! IN          CHAMP   K19: NOM DU CHAMP
! IN     ICHIN     I : NUMERO DE LA CHARGE ASSOCIEE A CHAMP
! OUT    ICHOUT    I : NUMERO DE LA CHARGE ASSOCIEE A CHAMP
!
! --------------------------------------------------------------------------------------------------
!
!  CONVENTIONS :
! --------------
! SI ACTION='E', ICHOUT EST INUTILISE  (IBID)
!   SI ICHIN > 0 :ICHIN EST BIEN LA NUMERO DE LA CHARGE ASSOCIEE A CHAMP
!                 SI IL EXISTE UNE FONC_MULT, ELLE LUI SERA APPLIQUEE.
!   SI ICHIN = -1 :CHAMP N'EST ASSOCIE A AUCUNE FONC_MULT
!                 (CHARGEMENT DE DILATATION PAR EXEMPLE)
!   SI ICHIN = -2 :CHAMP EST "BIDON" : NUL OU INEXISTANT
!   SI ICHIN = 0 OU ICHIN < -2  : ERREUR FATALE.
!
! SI ACTION='L', ICHIN EST INUTILISE  (IBID)
!    SI ON A FAIT AU PREALABLE UN LIEN AVEC CHAMP, ON REND ICHOUT
!    SINON ON REND ICHOUT=0
!    ICHOUT=0  VEUT DIRE QUE LE CHAMP N'EST PAS RENSEIGNE :
!         SOIT IL NE L'A JAMAIS ETE, SOIT IL A ETE EFFACE (ACTION:'S')
!
! SI ACTION='S', ICHIN ET ICHOUT SONT INUTILISES  (IBID)
!    SI ON A FAIT AU PREALABLE UN LIEN AVEC CHAMP : OK
!    SINON  : ERREUR FATALE.
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: reptJv = '&&CORICH.REPT'
    character(len=24), parameter :: nuchJv = '&&CORICH.NUCH'
    character(len=24), parameter :: reptmp = '&&CORICH.REPTMP'
    integer(kind=8), parameter :: nbFieldInit = 50
    integer(kind=8) :: iret, nbField, nbFieldMaxi, iField, fieldIndx
    character(len=24) :: field, fieldCopy
    integer(kind=8) :: ichin, ichout
    integer(kind=8), pointer :: nuch(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    ichin = 0
    ichout = 0
    field = fieldZ(1:19)

! - Create objects if don't exist
    call jeexin(reptJv, iret)
    if (iret .eq. 0) then
        call jecreo(reptJv, 'V N K24')
        call jeecra(reptJv, 'NOMMAX', nbFieldInit)
        call jedetr(nuchJv)
        call wkvect(nuchJv, 'V V I', nbFieldInit, vi=nuch)
    end if

! - Upgrade objects
    call jelira(reptJv, 'NOMMAX', nbFieldMaxi)
    call jelira(reptJv, 'NOMUTI', nbField)
    if (nbField .gt. nbFieldMaxi-1) then
        call juveca(nuchJv, 2*nbFieldMaxi)
        call jedupo(reptJv, 'V', reptmp, .false._1)
        call jedetr(reptJv)
        call jecreo(reptJv, 'V N K24')
        call jeecra(reptJv, 'NOMMAX', 2*nbFieldMaxi)
        do iField = 1, nbField
            call jenuno(jexnum(reptmp, iField), fieldCopy)
            call jecroc(jexnom(reptJv, fieldCopy))
        end do
        call jedetr(reptmp)
    end if

! - Create link (E)
    if (actionZ .eq. 'E') then
        ichin = ichin_
        ASSERT(ichin .ne. 0)
        ASSERT(ichin .ge. -2)
        call jenonu(jexnom(reptJv, field), fieldIndx)
        if (fieldIndx .eq. 0) then
            call jecroc(jexnom(reptJv, field))
        end if
        call jenonu(jexnom(reptJv, field), fieldIndx)
        call jeveuo(nuchJv, 'E', vi=nuch)
        nuch(fieldIndx) = ichin

! - Read link (E)
    else if (actionZ .eq. 'L') then
        call jenonu(jexnom(reptJv, field), fieldIndx)
        if (fieldIndx .eq. 0) then
            ichout = 0
        else
            call jeveuo(nuchJv, 'L', vi=nuch)
            ichout = nuch(fieldIndx)
        end if
        ichout_ = ichout

! - Suppress link (S)
    else if (actionZ .eq. 'S') then
        call jenonu(jexnom(reptJv, field), fieldIndx)
        ASSERT(fieldIndx .gt. 0)
        call jeveuo(nuchJv, 'E', vi=nuch)
        nuch(fieldIndx) = 0

    else
        ASSERT(ASTER_FALSE)

    end if
!
    call jedema()
end subroutine
