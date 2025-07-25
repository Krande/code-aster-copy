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

subroutine rsagsd(nomsd, ilong)
    implicit none
#include "jeveux.h"
#include "asterc/isnnem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: ilong
    character(len=*) :: nomsd
! person_in_charge: jacques.pellet at edf.fr
!       REDIMENSIONNEMENT D'UNE STRUCTURE DE DONNEES "RESULTAT-COMPOSE"
!       (LA TAILLE EST DOUBLEE SI LA LONGEUR VAUT 0)
!       LA SD RESTE INCHANGEE SI ELLE EXISTE ET SI LA TAILLE DEMANDEE
!       EST INFERIEURE OU EGALE A L'ACTUELLE
! ----------------------------------------------------------------------
! IN  : NOMSD  : NOM DE LA STRUCTURE "RESULTAT" A AGRANDIR
! IN  : ILONG  : NOUVELLE LONGUEUR DE LA S.D.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: iundef, iret, nbcham, nbordr, nborlu, newnb, neword, neworl
    integer(kind=8) :: jtachg, jordrg, i, j, k, jordrv, jpara
    real(kind=8) :: rundef
    integer(kind=8) :: n1, n2, kk, ier1
    character(len=24) :: nomobj
    character(len=19) :: nomd2
    character(len=24), pointer :: tach(:) => null()
! ----------------------------------------------------------------------
    call jemarq()
    nomd2 = nomsd
    rundef = r8vide()
    iundef = isnnem()
!
    if (ilong .lt. 0) then
        call utmess('F', 'UTILITAI4_29')
    end if
    call jeexin(nomd2//'.DESC', iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI_40', sk=nomd2)
    end if
!
    call jelira(nomd2//'.DESC', 'NOMMAX', nbcham)
    call jelira(nomd2//'.ORDR', 'LONMAX', nbordr)
    call jelira(nomd2//'.ORDR', 'LONUTI', nborlu)
    if (ilong .eq. 0) then
        newnb = 2*nbordr
    else
        newnb = ilong
    end if
    if (newnb .le. nbordr) goto 999
    neword = min(newnb, nbordr)
    neworl = min(newnb, nborlu)
!
!
!    -- LE .DESC, .NOVA, .TAVA ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     CES OBJETS NE SONT PAS MODIFIES
!
!
!
!     -- LE .TACH ET LE .ORDR ---
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    call jeveuo(nomd2//'.TACH', 'L', jtachg)
    AS_ALLOCATE(vk24=tach, size=neword*nbcham)
    call jeveuo(nomd2//'.ORDR', 'L', jordrg)
    call wkvect('&&RSAGSD.ORDR', 'V V I', max(neworl, 1), jordrv)
    do i = 0, neworl-1
        zi(jordrv+i) = zi(jordrg+i)
    end do
    do i = 0, nbcham-1
        do j = 0, neword-1
            tach(1+j+i*neword) = zk24(jtachg+j+i*nbordr)
        end do
    end do
    call jedetr(nomd2//'.TACH')
    call jedetr(nomd2//'.ORDR')
    call jecrec(nomd2//'.TACH', 'G V K24', 'NU', 'CONTIG', 'CONSTANT', &
                nbcham)
    call jeecra(nomd2//'.TACH', 'LONMAX', newnb)
    call jeveuo(nomd2//'.TACH', 'E', jtachg)
    do k = 1, nbcham
        call jecroc(jexnum(nomd2//'.TACH', k))
    end do
!
    call wkvect(nomd2//'.ORDR', 'G V I', newnb, jordrg)
    call jeecra(nomd2//'.ORDR', 'LONUTI', neworl)
!
    do i = 0, neworl-1
        zi(jordrg+i) = zi(jordrv+i)
    end do
    do i = 0, nbcham-1
        do j = 0, neword-1
            zk24(jtachg+j+i*newnb) = tach(1+j+i*neword)
        end do
    end do
!
    AS_DEALLOCATE(vk24=tach)
    call jedetr('&&RSAGSD.ORDR')
!
!
!
!     -- LES PARAMETRES :
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!
    nomobj = nomd2//'.RSPR'
    call jeexin(nomobj, ier1)
    if (ier1 .gt. 0) then
        call jelira(nomobj, 'LONMAX', n1)
        n2 = n1/nbordr
        ASSERT(n1 .eq. n2*nbordr)
        call juveca(nomobj, n2*newnb)
        call jeveuo(nomobj, 'E', jpara)
        do kk = n2*nbordr+1, n2*newnb
            zr(jpara-1+kk) = rundef
        end do
    end if
!
!
    nomobj = nomd2//'.RSPC'
    call jeexin(nomobj, ier1)
    if (ier1 .gt. 0) then
        call jelira(nomobj, 'LONMAX', n1)
        n2 = n1/nbordr
        call juveca(nomobj, n2*newnb)
        call jeveuo(nomobj, 'E', jpara)
        do kk = n2*nbordr+1, n2*newnb
            zc(jpara-1+kk) = dcmplx(rundef, rundef)
        end do
    end if
!
!
    nomobj = nomd2//'.RSPI'
    call jeexin(nomobj, ier1)
    if (ier1 .gt. 0) then
        call jelira(nomobj, 'LONMAX', n1)
        n2 = n1/nbordr
        call juveca(nomobj, n2*newnb)
        call jeveuo(nomobj, 'E', jpara)
        do kk = n2*nbordr+1, n2*newnb
            zi(jpara-1+kk) = iundef
        end do
    end if
!
!
    nomobj = nomd2//'.RSP8'
    call jeexin(nomobj, ier1)
    if (ier1 .gt. 0) then
        call jelira(nomobj, 'LONMAX', n1)
        n2 = n1/nbordr
        call juveca(nomobj, n2*newnb)
    end if
!
!
    nomobj = nomd2//'.RS16'
    call jeexin(nomobj, ier1)
    if (ier1 .gt. 0) then
        call jelira(nomobj, 'LONMAX', n1)
        n2 = n1/nbordr
        call juveca(nomobj, n2*newnb)
    end if
!
!
    nomobj = nomd2//'.RS24'
    call jeexin(nomobj, ier1)
    if (ier1 .gt. 0) then
        call jelira(nomobj, 'LONMAX', n1)
        n2 = n1/nbordr
        call juveca(nomobj, n2*newnb)
    end if
!
!
    nomobj = nomd2//'.RS32'
    call jeexin(nomobj, ier1)
    if (ier1 .gt. 0) then
        call jelira(nomobj, 'LONMAX', n1)
        n2 = n1/nbordr
        call juveca(nomobj, n2*newnb)
    end if
!
!
    nomobj = nomd2//'.RS80'
    call jeexin(nomobj, ier1)
    if (ier1 .gt. 0) then
        call jelira(nomobj, 'LONMAX', n1)
        n2 = n1/nbordr
        call juveca(nomobj, n2*newnb)
    end if
!
999 continue
!
    call jedema()
!
end subroutine
