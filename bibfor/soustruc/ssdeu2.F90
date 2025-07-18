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
subroutine ssdeu2(nval, iliste, nvalap)
    implicit none
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nval, iliste(nval), nvalap
! ----------------------------------------------------------------------
!     BUT:
!        - ENLEVER LES DOUBLONS D'UNE LISTE D'ENTIERS.
!          LES ENTIERS SONT RETASSES VERS LE DEBUT DE LA LISTE.
!        - ENLEVER LES "ZERO".
!
!     IN:
!        NVAL  :  NOMBRE D'ENTIERS A PRENDRE EN COMPTE.
!     IN/OUT:
!        ILISTE:  LISTE DES ENTIERS.
!                 EN SORTIE, ELLE EST RETASSEE. LA FIN DE LA LISTE EST
!                 MISE A ZERO.
!     OUT:
!        NVALAP:  NOMBRE D'ENTIERS DIFFERENTS TROUVES DANS LA LISTE.
!
! ----------------------------------------------------------------------
!
!
!     -- L'OBJET DE TRAVAIL "&&SSDEU2.WK1" CONTIENDRA DES "1" AU NIVEAU
!        DES ENTIERS A ELIMINER.
!     ---------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iawk1, idecal, iret, j, ndim
!-----------------------------------------------------------------------
    call jemarq()
    call jeexin('&&SSDEU2.WK1', iret)
    if (iret .eq. 0) then
        ndim = max(1000, 2*nval)
        call wkvect('&&SSDEU2.WK1', 'V V I', ndim, iawk1)
    else
        call jelira('&&SSDEU2.WK1', 'LONMAX', ndim)
        if (ndim .lt. nval) then
            call jedetr('&&SSDEU2.WK1')
            call wkvect('&&SSDEU2.WK1', 'V V I', 2*nval, iawk1)
        else
            call jeveuo('&&SSDEU2.WK1', 'E', iawk1)
        end if
    end if
!
!     -- MISE A ZERO DE "&&SSDEU2.WK1":
!     ---------------------------------
    do i = 1, nval
        zi(iawk1-1+i) = 0
    end do
!
!     -- MISE A "1" PARTIELLE DE  "&&SSDEU2.WK1":
!     -------------------------------------------
    nvalap = nval
    do 1, i = 1, nval
    if (iliste(i) .eq. 0) then
        zi(iawk1-1+i) = 1
        nvalap = nvalap-1
        goto 1
    end if
    do 2, j = 1, i-1
    if (iliste(j) .eq. iliste(i)) then
        zi(iawk1-1+i) = 1
        nvalap = nvalap-1
        goto 1
    end if
2   continue
1   end do
!
!
!     -- RETASSAGE DE LA LISTE:
!     -------------------------
    idecal = 0
    do 3, i = 1, nval
    if (zi(iawk1-1+i) .eq. 1) then
        idecal = idecal+1
    else
        iliste(i-idecal) = iliste(i)
    end if
3   end do
!
!     -- ON COMPLETE PAR DES ZERO:
!     ----------------------------
    do 4, i = nvalap+1, nval
        iliste(i) = 0
4   end do
!
!
    call jedema()
    end subroutine
