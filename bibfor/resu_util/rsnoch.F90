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

subroutine rsnoch(nomsd, nomsy, iordr)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutrg.h"
#include "asterfort/sdmpic.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: iordr
    character(len=*) :: nomsd, nomsy
! person_in_charge: jacques.pellet at edf.fr
!
!  BUT : "NOTER" UN CHAMP DANS UNE SD_RESULTAT
!        ON VERIFIE QUE :
!           - LA PLACE EST LICITE (NOMSY OK ET IORDR<=NBORDR_MAX)
!           - LE CHAMP QUI VA ETRE "NOTE" DANS NOMSD
!             (REGLE DE NOMMAGE DE RSUTCH.F)
!             EXISTE REELLEMENT (EXISD.F)
! ----------------------------------------------------------------------
! IN  : NOMSD  : NOM DE LA STRUCTURE "RESULTAT"
! IN  : NOMSY  : NOM SYMBOLIQUE DU CHAMP A NOTER.
! IN  : IORDR  : NUMERO D'ORDRE DU CHAMP A NOTER.
! ----------------------------------------------------------------------
!
    character(len=16) :: noms2
    character(len=19) :: nomd2, chnote
    character(len=24) :: valk(2)
    character(len=8) :: repk
    integer(kind=8) :: normax, iretou, nordr, irang, iret, ibid, jtach
    integer(kind=8), pointer :: ordr(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    noms2 = nomsy
    nomd2 = nomsd
!
!
!     -- CALCUL ET VALIDATION DU NUMERO DE RANGEMENT :IRANG
!     -----------------------------------------------------
    call jelira(nomd2//'.ORDR', 'LONMAX', normax)
    call rsutrg(nomd2, iordr, iretou, nordr)
    if (iretou .eq. 0) then
        irang = nordr+1
        if (irang .gt. normax) then
            call utmess('F', 'UTILITAI4_42')
        end if
        call jeecra(nomd2//'.ORDR', 'LONUTI', irang)
        call jeveuo(nomd2//'.ORDR', 'E', vi=ordr)
!       -- ON VERIFIE QUE LE NOUVEAU IORDR EST SUPERIEUR
!          AU DERNIER IORDR DEJA STOCKE (IORDR CROISSANTS) :
        if (irang .gt. 1) then
            if (ordr(irang-1) >= iordr) then
                call utmess('F', 'UTILITAI6_81', sk=nomsd, ni=2, vali=[iordr, nordr])
            end if
        end if
        ordr(irang) = iordr
    else
        irang = iretou
    end if
!
!
!     -- ON VERIFIE LE NOM SYMBOLIQUE :
!     -------------------------------------------
    call jenonu(jexnom(nomd2//'.DESC', noms2), iret)
    if (iret .eq. 0) then
        valk(1) = noms2
        valk(2) = nomd2
        call utmess('F', 'UTILITAI4_43', nk=2, valk=valk)
    end if
!
!
!     -- CHNOTE : NOM QUE DOIT AVOIR LE CHAMP A NOTER :
!        (REGLE DE NOMMAGE DE RSUTCH.F)
!     -------------------------------------------------
    call rsexch(' ', nomd2, noms2, iordr, chnote, &
                iret)
!
!     -- ON VERIFIE L'EXISTENCE DE CHNOTE :
!     -------------------------------------------
    if (iret .eq. 100) then
        valk(1) = noms2
        valk(2) = chnote
        call utmess('F', 'UTILITAI_50', nk=2, valk=valk)
    end if
    ASSERT(iret .eq. 0)
!
!
!     --- ON STOCKE LE NOM DU CHAMP :
!     ------------------------------
    call jenonu(jexnom(nomd2//'.DESC', noms2), ibid)
    call jeveuo(jexnum(nomd2//'.TACH', ibid), 'E', jtach)
!
    zk24(jtach+irang-1) (1:19) = chnote
!
!
!     -- SI LE CHAMP EST UN CHAM_ELEM MPI_INCOMPLET, ON LE COMPLETE:
!     --------------------------------------------------------------
    call dismoi('TYPE_CHAMP', chnote, 'CHAMP', repk=repk)
    if (repk(1:2) .eq. 'EL') call sdmpic('CHAM_ELEM', chnote)
!
    call jedema()
end subroutine
