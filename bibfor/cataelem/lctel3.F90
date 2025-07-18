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

subroutine lctel3()
    implicit none
!
! ----------------------------------------------------------------------
!     BUT: CREER APOSTERIORI L'OBJET &CATA.TE.DIM_GEOM
!          QUI CONTIENT LA DIMENSION GEOMETRIQUE DES TYPE_ELEM
!          0 : LE TYPE_ELEM N'UTILISE PAS LA GRANDEUR "GEOM_R"
!          1 : LE TYPE_ELEM UTILISE LA CMP "X"
!          2 : LE TYPE_ELEM UTILISE LA CMP "Y"
!          3 : LE TYPE_ELEM UTILISE LA CMP "Z"
!
! ----------------------------------------------------------------------
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: dg
    character(len=16) :: nomte
    character(len=24) :: nomolo
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iadige, iamolo, icode, igd, igdgeo, iml
    integer(kind=8) :: inocmp, ite, ix, iy, iz, k, nbcmp
    integer(kind=8) :: nbdg, nbec, nbml, nbpt, nbte
!-----------------------------------------------------------------------
    call jemarq()
    call jelira('&CATA.TE.NOMTE', 'NOMMAX', nbte)
    call wkvect('&CATA.TE.DIM_GEOM', 'G V I', nbte, iadige)
    call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), igdgeo)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', 'GEOM_R'), 'L', inocmp)
    call jelira(jexnom('&CATA.GD.NOMCMP', 'GEOM_R'), 'LONMAX', nbcmp)
    ix = indik8(zk8(inocmp), 'X', 1, nbcmp)
    iy = indik8(zk8(inocmp), 'Y', 1, nbcmp)
    iz = indik8(zk8(inocmp), 'Z', 1, nbcmp)
    call dismoi('NB_EC', 'GEOM_R', 'GRANDEUR', repi=nbec)
    ASSERT(nbec .le. 1)
!
!     - BOUCLE SUR TOUS LES MODES LOCAUX DES CATALOGUES :
    call jelira('&CATA.TE.NOMMOLOC', 'NOMMAX', nbml)
    do iml = 1, nbml
        call jeveuo(jexnum('&CATA.TE.MODELOC', iml), 'L', iamolo)
        icode = zi(iamolo-1+1)
        igd = zi(iamolo-1+2)
        if (igd .ne. igdgeo) goto 1
        if (icode .gt. 3) goto 1
!
        call jenuno(jexnum('&CATA.TE.NOMMOLOC', iml), nomolo)
        nomte = nomolo(1:16)
        call jenonu(jexnom('&CATA.TE.NOMTE', nomte), ite)
!
        nbpt = zi(iamolo-1+4)
        if (nbpt .ge. 10000) then
            nbdg = nbpt-10000
        else
            nbdg = 1
        end if
!
        do k = 1, nbdg
            dg = zi(iamolo-1+4+k)
            if (exisdg([dg], ix)) zi(iadige-1+ite) = max(1, zi(iadige-1+ite))
            if (exisdg([dg], iy)) zi(iadige-1+ite) = max(2, zi(iadige-1+ite))
            if (exisdg([dg], iz)) zi(iadige-1+ite) = max(3, zi(iadige-1+ite))
        end do
1       continue
    end do
!
!
    call jedema()
end subroutine
