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
subroutine ssdeu1(motcle, noma, nbno, iliste)
    implicit none
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/getvem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    character(len=*) :: motcle
    character(len=8) :: noma
    integer(kind=8) :: nbno, iliste(*)
! ----------------------------------------------------------------------
!     BUT:
!        - TRAITER LES MOTS CLEFS "GROUP_NO" ET "NOEUD" DE
!          EXTERIEUR DE LA COMMANDE MACR_ELEM_STAT.
!        - COMPTER LES NOEUDS TROUVES , EN ETABLIR LA LISTE.
!
!     IN:
!        MOTCLE: 'NOMBRE' --> ON REND LE NOMBRE DE NOEUDS UNIQUEMENT.
!                'LISTE ' --> ON REND LE NOMBRE ET LA LISTE.
!        NOMA  : NOM DU MAILLAGE
!     OUT:
!        NBNO  :  NOMBRE DE NOEUDS TROUVES.
!        ILISTE:  LISTE DES NUMEROS DE NOEUDS TROUVES
!                 (SI MOTCLE='LISTE' SEULEMENT.)
!
    character(len=8) :: kbi81
    character(len=24) :: valk(2)
! ----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iagpno, iawk1, ibid, ico, ii, iret
    integer(kind=8) :: n1, n2, n3, n4, ndim
!-----------------------------------------------------------------------
    call jemarq()
    nbno = 0
!
!
!
    call jeexin('&&SSDEU1.WK1', iret)
    if (iret .le. 0) then
        ndim = 100
        call wkvect('&&SSDEU1.WK1', 'V V K8', ndim, iawk1)
    else
        call jelira('&&SSDEU1.WK1', 'LONMAX', ndim)
        call jeveuo('&&SSDEU1.WK1', 'E', iawk1)
    end if
!
!
!     --CAS NOEUD:
!     ------------
    call getvem(noma, 'NOEUD', 'EXTERIEUR', 'NOEUD', 1, &
                0, kbi81, n1)
    if (n1 .ne. 0) then
        n3 = -n1
        if (ndim .lt. n3) then
            call jedetr('&&SSDEU1.WK1')
            call wkvect('&&SSDEU1.WK1', 'V V K8', 2*n3, iawk1)
        end if
        call getvem(noma, 'NOEUD', 'EXTERIEUR', 'NOEUD', 1, &
                    n3, zk8(iawk1), ibid)
        nbno = nbno+n3
        if (motcle .eq. 'LISTE') then
            do i = 1, n3
                iliste(i) = char8_to_int(zk8(iawk1-1+i))
                if (iliste(i) .eq. 0) then
                    valk(1) = zk8(iawk1-1+i)
                    valk(2) = noma
                    call utmess('F', 'SOUSTRUC_48', nk=2, valk=valk)
                end if
            end do
        end if
    end if
!
!
!     --CAS GROUP_NO:
!     ---------------
    call getvem(noma, 'GROUP_NO', 'EXTERIEUR', 'GROUP_NO', 1, &
                0, kbi81, n2)
    if (n2 .ne. 0) then
        n3 = -n2
        call jedetr('&&SSDEU1.WK1')
        call wkvect('&&SSDEU1.WK1', 'V V K24', 2*n3, iawk1)
        call getvem(noma, 'GROUP_NO', 'EXTERIEUR', 'GROUP_NO', 1, &
                    n3, zk24(iawk1), ibid)
        ico = nbno
        do i = 1, n3
            call jeexin(jexnom(noma//'.GROUPENO', zk24(iawk1-1+i)), iret)
            if (iret .eq. 0) then
                valk(1) = zk24(iawk1-1+i)
                valk(2) = noma
                call utmess('F', 'SOUSTRUC_49', nk=2, valk=valk)
            end if
            call jelira(jexnom(noma//'.GROUPENO', zk24(iawk1-1+i)), 'LONMAX', n4)
            nbno = nbno+n4
            if (motcle .eq. 'LISTE') then
                call jeveuo(jexnom(noma//'.GROUPENO', zk24(iawk1-1+i)), 'L', iagpno)
                do ii = 1, n4
                    ico = ico+1
                    iliste(ico) = zi(iagpno-1+ii)
                end do
            end if
        end do
    end if
!
!
!
    call jedema()
end subroutine
