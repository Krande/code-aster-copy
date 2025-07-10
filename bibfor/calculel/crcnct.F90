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

subroutine crcnct(base, nomch, mailla, gd, nbcmp, &
                  licmp, rcmp)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/afchno.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
!
    integer(kind=8) :: nbcmp
    character(len=*) :: base, nomch, mailla, gd
    character(len=*) :: licmp(nbcmp)
    real(kind=8) :: rcmp(nbcmp)
! ----------------------------------------------------------------------
!     BUT:
!     ----
!      CREER UN CHAM_NO A REPRESENTATION CONSTANTE D'UNE GRANDEUR
!      DONNEE AVEC DES COMPOSANTES IDENTIQUES SUR TOUS LES NOEUDS.
!
!     ENTREES:
!     --------
!     BASE  (K1)  : BASE OU L'ON VEUT CREER LE CHAM_NO
!     MAILLA(K8)  : NOM DU MAILLAGE SUR LEQUEL ON VEUT CREER LE CHAM_NO
!     GD    (K8)  : NOM DE LA GRANDEUR ASSICIEE AU CHAMP
!     NBCMP (I)   : NOMBRE DE CMPS DE LICMP(*)
!     LICMP (LK8) : LISTE DES CMPS QUE L'ON VEUT SUR TOUS LES NOEUDS
!     RCMP  (LR)  : VALEURS DES CMPS QUE L'ON VEUT EN TOUS LES NOEUDS
!
!     SORTIES:
!     --------
!      LE CHAM_NO NOMCH EST CREE.
!
! ----------------------------------------------------------------------
    character(len=19) :: ch19
    character(len=8) :: maill2, gd2, nocmp
    character(len=1) :: bas2
    character(len=4) :: tysca
    character(len=24) :: valk(2)
!
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iancmp, iavale, icmp, iec
    integer(kind=8) :: igd, ino, itrou, nbcmp2, nbno
    integer(kind=8) :: nec, jj
    integer(kind=8), pointer :: desc(:) => null()
    integer(kind=8), pointer :: nbca(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    gd2 = gd
    maill2 = mailla
    ch19 = nomch
    bas2 = base

!
!
!     VERIFICATION DES ARGUMENTS D'APPEL :
!     ------------------------------------
    call jenonu(jexnom('&CATA.GD.NOMGD', gd2), igd)
    if (igd .eq. 0) then
        call utmess('F', 'CALCULEL2_21', sk=gd2)
    end if
    call jeveuo(jexnum('&CATA.GD.NOMCMP', igd), 'L', iancmp)
    call jelira(jexnum('&CATA.GD.NOMCMP', igd), 'LONMAX', nbcmp2)
    call dismoi('NB_EC', gd2, 'GRANDEUR', repi=nec)
    call dismoi('TYPE_SCA', gd2, 'GRANDEUR', repk=tysca)
    if (tysca(1:1) .ne. 'R') then
        call utmess('F', 'CALCULEL2_23', sk=gd2)
    end if
!
!
!     ALLOCATION DU CHAM_NO :
!     -----------------------
    call dismoi('NB_NO_MAILLA', maill2, 'MAILLAGE', repi=nbno)
    call wkvect('&&CRCNCT.NCMPMX_AFFE', 'V V I ', nbno, vi=nbca)
    call wkvect('&&CRCNCT.DESC_NOEUD', 'V V I', nec*nbno, vi=desc)

    do ino = 1, nbno
        nbca(ino) = nbcmp
        do icmp = 1, nbcmp
            nocmp = licmp(icmp)
            itrou = indik8(zk8(iancmp), nocmp, 1, nbcmp2)
            if (itrou .ne. 0) then
                iec = (itrou-1)/30+1
                jj = itrou-30*(iec-1)
                desc((ino-1)*nec+iec) = ior(desc((ino-1)*nec+iec), 2**jj)
            else
                valk(1) = nocmp
                valk(2) = gd2
                call utmess('F', 'CALCULEL2_22', nk=2, valk=valk)
            end if
        end do
    end do

    call afchno(ch19, bas2, gd, maill2, nbno, nbca, desc, nbcmp*nbno, 'R')
    call jeveuo(ch19//'.VALE', 'E', iavale)
!
!     OBJET: .VALE
!     ------------
    do icmp = 1, nbcmp
        do ino = 1, nbno
            zr(iavale-1+(ino-1)*nbcmp+icmp) = rcmp(icmp)
        end do
    end do
!
    call jedetr('&&CRCNCT.NCMPMX_AFFE')
    call jedetr('&&CRCNCT.DESC_NOEUD')
!
    call jedema()
end subroutine
