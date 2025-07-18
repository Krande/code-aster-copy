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

subroutine giinco()
    implicit none
!
! ----------------------------------------------------------------------
!     BUT: CREER LA COLLECTION '&&GILIRE.CORR_GIBI_ASTER'
!          QUI DONNE LA PERMUTATION DES NOEUDS DES MAILLES
!          (GIBI--> ASTER)
!          CETTE COLLECTION EST UTILISEE DANS GIECMA.
!
! ----------------------------------------------------------------------
!
!     VARIABLES LOCALES:
#include "jeveux.h"
!
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
    character(len=24) :: nomcol
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iacorr
!-----------------------------------------------------------------------
    call jemarq()
    nomcol = '&&GILIRE.CORR_GIBI_ASTER'
!
    call jecrec(nomcol, 'V V I', 'NO', 'CONTIG', 'VARIABLE', &
                17)
!
!     ON DIMENSIONNE LA COLLECTION:
!     -----------------------------
    call jecroc(jexnom(nomcol, 'POI1'))
    call jeecra(jexnom(nomcol, 'POI1'), 'LONMAX', 1)
    call jecroc(jexnom(nomcol, 'SEG2'))
    call jeecra(jexnom(nomcol, 'SEG2'), 'LONMAX', 2)
    call jecroc(jexnom(nomcol, 'SEG3'))
    call jeecra(jexnom(nomcol, 'SEG3'), 'LONMAX', 3)
    call jecroc(jexnom(nomcol, 'TRI3'))
    call jeecra(jexnom(nomcol, 'TRI3'), 'LONMAX', 3)
    call jecroc(jexnom(nomcol, 'TRI6'))
    call jeecra(jexnom(nomcol, 'TRI6'), 'LONMAX', 6)
    call jecroc(jexnom(nomcol, 'QUA4'))
    call jeecra(jexnom(nomcol, 'QUA4'), 'LONMAX', 4)
    call jecroc(jexnom(nomcol, 'QUA8'))
    call jeecra(jexnom(nomcol, 'QUA8'), 'LONMAX', 8)
    call jecroc(jexnom(nomcol, 'QUA9'))
    call jeecra(jexnom(nomcol, 'QUA9'), 'LONMAX', 9)
    call jecroc(jexnom(nomcol, 'CUB8'))
    call jeecra(jexnom(nomcol, 'CUB8'), 'LONMAX', 8)
    call jecroc(jexnom(nomcol, 'CU20'))
    call jeecra(jexnom(nomcol, 'CU20'), 'LONMAX', 20)
    call jecroc(jexnom(nomcol, 'CU27'))
    call jeecra(jexnom(nomcol, 'CU27'), 'LONMAX', 27)
    call jecroc(jexnom(nomcol, 'PRI6'))
    call jeecra(jexnom(nomcol, 'PRI6'), 'LONMAX', 6)
    call jecroc(jexnom(nomcol, 'PR15'))
    call jeecra(jexnom(nomcol, 'PR15'), 'LONMAX', 16)
    call jecroc(jexnom(nomcol, 'TET4'))
    call jeecra(jexnom(nomcol, 'TET4'), 'LONMAX', 4)
    call jecroc(jexnom(nomcol, 'TE10'))
    call jeecra(jexnom(nomcol, 'TE10'), 'LONMAX', 10)
    call jecroc(jexnom(nomcol, 'PYR5'))
    call jeecra(jexnom(nomcol, 'PYR5'), 'LONMAX', 5)
    call jecroc(jexnom(nomcol, 'PY13'))
    call jeecra(jexnom(nomcol, 'PY13'), 'LONMAX', 13)
!
!     POI1:
!     -----
    call jeveuo(jexnom(nomcol, 'POI1'), 'E', iacorr)
    zi(iacorr-1+1) = 1
!
!     SEG2:
!     -----
    call jeveuo(jexnom(nomcol, 'SEG2'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 2
!
!     SEG3:
!     -----
    call jeveuo(jexnom(nomcol, 'SEG3'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 3
    zi(iacorr-1+3) = 2
!
!     TRI3:
!     -----
    call jeveuo(jexnom(nomcol, 'TRI3'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 2
    zi(iacorr-1+3) = 3
!
!     TRI6:
!     -----
    call jeveuo(jexnom(nomcol, 'TRI6'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 3
    zi(iacorr-1+3) = 5
    zi(iacorr-1+4) = 2
    zi(iacorr-1+5) = 4
    zi(iacorr-1+6) = 6
!
!     QUA4:
!     -----
    call jeveuo(jexnom(nomcol, 'QUA4'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 2
    zi(iacorr-1+3) = 3
    zi(iacorr-1+4) = 4
!
!     QUA8:
!     -----
    call jeveuo(jexnom(nomcol, 'QUA8'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 3
    zi(iacorr-1+3) = 5
    zi(iacorr-1+4) = 7
    zi(iacorr-1+5) = 2
    zi(iacorr-1+6) = 4
    zi(iacorr-1+7) = 6
    zi(iacorr-1+8) = 8
!
!     QUA9:
!     -----
    call jeveuo(jexnom(nomcol, 'QUA9'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 3
    zi(iacorr-1+3) = 5
    zi(iacorr-1+4) = 7
    zi(iacorr-1+5) = 2
    zi(iacorr-1+6) = 4
    zi(iacorr-1+7) = 6
    zi(iacorr-1+8) = 8
    zi(iacorr-1+9) = 9
!
!     CUB8:
!     -----
    call jeveuo(jexnom(nomcol, 'CUB8'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 2
    zi(iacorr-1+3) = 3
    zi(iacorr-1+4) = 4
    zi(iacorr-1+5) = 5
    zi(iacorr-1+6) = 6
    zi(iacorr-1+7) = 7
    zi(iacorr-1+8) = 8
!
!     CU20:
!     -----
    call jeveuo(jexnom(nomcol, 'CU20'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 3
    zi(iacorr-1+3) = 5
    zi(iacorr-1+4) = 7
    zi(iacorr-1+5) = 13
    zi(iacorr-1+6) = 15
    zi(iacorr-1+7) = 17
    zi(iacorr-1+8) = 19
    zi(iacorr-1+9) = 2
    zi(iacorr-1+10) = 4
    zi(iacorr-1+11) = 6
    zi(iacorr-1+12) = 8
    zi(iacorr-1+13) = 9
    zi(iacorr-1+14) = 10
    zi(iacorr-1+15) = 11
    zi(iacorr-1+16) = 12
    zi(iacorr-1+17) = 14
    zi(iacorr-1+18) = 16
    zi(iacorr-1+19) = 18
    zi(iacorr-1+20) = 20
!
!     CU27:
!     -----
    call jeveuo(jexnom(nomcol, 'CU27'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 3
    zi(iacorr-1+3) = 5
    zi(iacorr-1+4) = 7
    zi(iacorr-1+5) = 13
    zi(iacorr-1+6) = 15
    zi(iacorr-1+7) = 17
    zi(iacorr-1+8) = 19
    zi(iacorr-1+9) = 2
    zi(iacorr-1+10) = 4
    zi(iacorr-1+11) = 6
    zi(iacorr-1+12) = 8
    zi(iacorr-1+13) = 9
    zi(iacorr-1+14) = 10
    zi(iacorr-1+15) = 11
    zi(iacorr-1+16) = 12
    zi(iacorr-1+17) = 14
    zi(iacorr-1+18) = 16
    zi(iacorr-1+19) = 18
    zi(iacorr-1+20) = 20
    zi(iacorr-1+21) = 25
    zi(iacorr-1+22) = 21
    zi(iacorr-1+23) = 22
    zi(iacorr-1+24) = 23
    zi(iacorr-1+25) = 24
    zi(iacorr-1+26) = 26
    zi(iacorr-1+27) = 27
!
!     PRI6:
!     -----
    call jeveuo(jexnom(nomcol, 'PRI6'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 2
    zi(iacorr-1+3) = 3
    zi(iacorr-1+4) = 4
    zi(iacorr-1+5) = 5
    zi(iacorr-1+6) = 6
!
!     PR15:
!     -----
    call jeveuo(jexnom(nomcol, 'PR15'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 3
    zi(iacorr-1+3) = 5
    zi(iacorr-1+4) = 10
    zi(iacorr-1+5) = 12
    zi(iacorr-1+6) = 14
    zi(iacorr-1+7) = 2
    zi(iacorr-1+8) = 4
    zi(iacorr-1+9) = 6
    zi(iacorr-1+10) = 7
    zi(iacorr-1+11) = 8
    zi(iacorr-1+12) = 9
    zi(iacorr-1+13) = 11
    zi(iacorr-1+14) = 13
    zi(iacorr-1+15) = 15
!
!     TET4:
!     -----
    call jeveuo(jexnom(nomcol, 'TET4'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 2
    zi(iacorr-1+3) = 3
    zi(iacorr-1+4) = 4
!
!     TE10:
!     -----
    call jeveuo(jexnom(nomcol, 'TE10'), 'E', iacorr)
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 3
    zi(iacorr-1+3) = 5
    zi(iacorr-1+4) = 10
    zi(iacorr-1+5) = 2
    zi(iacorr-1+6) = 4
    zi(iacorr-1+7) = 6
    zi(iacorr-1+8) = 7
    zi(iacorr-1+9) = 8
    zi(iacorr-1+10) = 9
!
!     PYR5:
!     -----
    call jeveuo(jexnom(nomcol, 'PYR5'), 'E', iacorr)
!     ON NE CONNAIT PAS ENCORE LA NUMEROTATION ASTER ...
    zi(iacorr-1+1) = 1
    zi(iacorr-1+2) = 2
    zi(iacorr-1+3) = 3
    zi(iacorr-1+4) = 4
    zi(iacorr-1+5) = 5
!
!     PY13:
!     -----
    call jeveuo(jexnom(nomcol, 'PY13'), 'E', iacorr)
!     ON NE CONNAIT PAS ENCORE LA NUMEROTATION ASTER ...
    zi(iacorr-1+1) = 1
    zi(iacorr-1+6) = 2
    zi(iacorr-1+2) = 3
    zi(iacorr-1+7) = 4
    zi(iacorr-1+3) = 5
    zi(iacorr-1+8) = 6
    zi(iacorr-1+4) = 7
    zi(iacorr-1+9) = 8
    zi(iacorr-1+10) = 9
    zi(iacorr-1+11) = 10
    zi(iacorr-1+12) = 11
    zi(iacorr-1+13) = 12
    zi(iacorr-1+5) = 13
!
!
!
    call jedema()
end subroutine
