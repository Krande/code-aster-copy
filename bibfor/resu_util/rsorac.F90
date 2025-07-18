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
subroutine rsorac(nomsd, acces, ival, rval, kval, &
                  cval, epsi, crit, nutrou, ndim, &
                  nbtrou)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxliis.h"
#include "asterfort/rsindi.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(out) :: nbtrou, nutrou(*)
    integer(kind=8), intent(in) :: ival, ndim
    real(kind=8), intent(in) :: rval, epsi
    character(len=*), intent(in) :: nomsd, acces, kval, crit
    complex(kind=8), intent(in) :: cval
! person_in_charge: jacques.pellet at edf.fr
!      RECUPERATION DU NUMERO D'ORDRE
!      D'UNE STRUCTURE DE DONNEES "SD_RESULTAT".
!      A PARTIR D'UNE VARIABLE D'ACCES.
!      ( CETTE ROUTINE FONCTIONNE AUSSI AVEC UN CONCEPT CHAMP_GD SI
!        RVAL,IVAL,..= 0 OU ACCES='DERNIER' ALORS NUTROU=1 )
! ----------------------------------------------------------------------
! IN  : NOMSD  : NOM DE LA SD "SD_RESULTAT".
!
! IN  : ACCES  : NOM SYMBOLIQUE DE LA VARIABLE D'ACCES.
!      OU BIEN : 'LONUTI','LONMAX','PREMIER','DERNIER','TOUT_ORDRE'
! ATTENTION : ACCES='LONUTI' OU 'LONMAX' NE RENVOIENT PAS UNE LISTE
!    DE NUMEROS D'ORDRE MAIS LEUR NOMBRE.
!    NBTROU=1 ET NUTROU= NOMBRE TROUVE !!
!
! IN  : IVAL   : VALEUR DE LA VARIABLE D'ACCES (SI ENTIER).
! IN  : RVAL   : VALEUR DE LA VARIABLE D'ACCES (SI REEL).
! IN  : CVAL   : VALEUR DE LA VARIABLE D'ACCES (SI COMPLEXE).
! IN  : EPSI   : PRECISION DE LA VARIABLE D'ACCES (RELATIVE/ABSOLUE).
! IN  : CRIT   : CRITERE DE PRECISION : 'RELATIF' OU 'ABSOLU'
!                (UNE VARIABLE D'ACCES EST DECLAREE VALIDE SI ELLE
!                                  SATISFAIT LE TEST RELATIF OU ABSOLU)
! IN  : NDIM   : DIMENSION DE LA LISTE NUTROU.
! OUT : NUTROU : LISTE DES NUMEROS D'ORDRE TROUVES.
! OUT : NBTROU : NOMBRE DE NUMEROS D'ORDRE TROUVES.
!              (SI LE NOMBRE TROUVE EST > NDIM, ON REND NBTROU=-NBTROU)
! ----------------------------------------------------------------------
    character(len=4) :: tysd, type, tysca
    character(len=8) :: nomobj, k8debu, k8maxi, k8ent
    character(len=16) :: acce2
    character(len=19) :: noms2
! ----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacces, iaobj, iatava, idebu, ier1
    integer(kind=8) :: ier2, iloty, imaxi, jordr, nbordr, nordr, numed
!
!-----------------------------------------------------------------------
    call jemarq()
    noms2 = nomsd
    acce2 = acces
!
!
!     --- CONCEPT CHAMP-GD
!     ----------------------------
    call jelira(noms2//'.DESC', 'DOCU', cval=tysd)
    if ((tysd .eq. 'CHNO') .or. (tysd .eq. 'CHML') .or. (tysd .eq. 'CART')) then
        if ((acce2 .eq. 'LONUTI') .or. (ival .eq. 0) .or. (rval .eq. 0.d0) .or. &
            (cval .eq. (0.d0, 0.d0))) then
            if (ndim .gt. 0) then
                nbtrou = 1
                nutrou(1) = 1
            else
                nbtrou = -1
            end if
        else
            call utmess('F', 'UTILITAI4_46')
        end if
        goto 20
    end if
!
!
!     --- CONCEPT RESULTAT
!     ----------------------------
    if (acce2 .eq. 'LONUTI') then
        if (ndim .gt. 0) then
            nbtrou = 1
            call jelira(noms2//'.ORDR', 'LONUTI', nutrou(1))
        else
            nbtrou = -1
        end if
        goto 20
!
    else if (acce2 .eq. 'LONMAX') then
        if (ndim .gt. 0) then
            nbtrou = 1
            call jelira(noms2//'.ORDR', 'LONMAX', nutrou(1))
        else
            nbtrou = -1
        end if
        goto 20
!
    else if (acce2 .eq. 'DERNIER') then
        if (ndim .gt. 0) then
            nbtrou = 1
            call jelira(noms2//'.ORDR', 'LONUTI', numed)
            if (numed .eq. 0) then
                nbtrou = 0
            else
                call jeveuo(noms2//'.ORDR', 'L', jordr)
                nutrou(1) = zi(jordr+numed-1)
            end if
        else
            nbtrou = -1
        end if
        goto 20
!
    else if (acce2 .eq. 'PREMIER') then
        if (ndim .gt. 0) then
            nbtrou = 1
            call jeveuo(noms2//'.ORDR', 'L', jordr)
            nutrou(1) = zi(jordr-1+1)
        else
            nbtrou = -1
        end if
        goto 20
!
    else if (acce2 .eq. 'TOUT_ORDRE') then
        call jelira(noms2//'.ORDR', 'LONUTI', nordr)
        if (nordr .le. ndim) then
            nbtrou = nordr
            call jeveuo(noms2//'.ORDR', 'L', jordr)
            do i = 1, nordr
                nutrou(i) = zi(jordr+i-1)
            end do
        else
            nbtrou = -nordr
            call jeveuo(noms2//'.ORDR', 'L', jordr)
            do i = 1, ndim
                nutrou(i) = zi(jordr+i-1)
            end do
        end if
        goto 20
!
    end if
!
    call jenonu(jexnom(noms2//'.NOVA', acce2), iacces)
    if (iacces .eq. 0) then
        call utmess('F', 'UTILITAI4_47', sk=acce2)
    end if
!
    call jeveuo(jexnum(noms2//'.TAVA', iacces), 'L', iatava)
    nomobj = zk8(iatava-1+1)
    k8maxi = zk8(iatava-1+3)
    call lxliis(k8maxi, imaxi, ier2)
    k8debu = zk8(iatava-1+2)
    call lxliis(k8debu, idebu, ier1)
    ASSERT(imaxi .gt. 0)
    ASSERT((idebu .gt. 0) .and. (idebu .le. imaxi))
!
    call jeveuo(noms2//'.ORDR', 'L', jordr)
    call jeveuo(noms2//nomobj, 'L', iaobj)
    call jelira(noms2//nomobj, 'TYPE', cval=type)
    call jelira(noms2//nomobj, 'LTYP', iloty)
    call jelira(noms2//'.ORDR', 'LONUTI', nbordr)
    call codent(iloty, 'G', k8ent)
    tysca = type(1:1)//k8ent(1:3)
!
    call rsindi(tysca, iaobj-1+idebu, imaxi, jordr, ival, &
                rval, kval, cval, epsi, crit, &
                nbordr, nbtrou, nutrou, ndim)
!
20  continue
    call jedema()
end subroutine
