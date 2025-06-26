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
subroutine detgnm(ma)
    implicit none
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cpclma.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: ma
!
!     BUT: - TRAITE LES MOTS CLES FACTEURS DETR_GROUP_MA ET
!            DETR_GROUP_NO DE L'OPERATEUR DEFI_GROUP.
!          - PERMET DE DETRUIRE DES GROUPES (NOEUDS OU MAILLES)
!
!     IN : MA : NOM DU MAILLAGE
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: n1, iocc, maxval, nbval, iret, nbtgp, numgm, nbgmde
    integer(kind=8) :: nbgp, nbmagp, jgm, i, j, j1, j2, ngp, ig
    parameter(ngp=2)
    character(len=8) :: k8b
    character(len=16) :: detr(2), group(2), ptrn(2)
    character(len=24) :: grp, gpptnm, nomgp
    character(len=24), pointer :: group_detr(:) => null()
    data detr/'DETR_GROUP_MA', 'DETR_GROUP_NO'/
    data group/'.GROUPEMA', '.GROUPENO'/
    data ptrn/'.PTRNOMMAI', '.PTRNOMNOE'/
!
    call jemarq()
!
    do ig = 1, ngp
        call getfac(detr(ig), n1)
        grp = '&&DETGNM'//group(ig)
        gpptnm = '&&DETGNM'//ptrn(ig)
        if (n1 .ne. 0) then
            call jeexin(ma//group(ig), iret)
            if (iret .eq. 0) goto 100
            call jelira(ma//group(ig), 'NUTIOC', nbtgp)
            call wkvect('&&DETGNM.GROUP', 'V V I', nbtgp, jgm)
            do i = 1, nbtgp
                zi(jgm+i-1) = 0
            end do
            do iocc = 1, n1
                maxval = 0
                call getvtx(detr(ig), 'NOM', iocc=iocc, nbval=maxval, vect=k8b, &
                            nbret=nbval)
                nbval = -nbval
                AS_ALLOCATE(vk24=group_detr, size=nbval)
                call getvtx(detr(ig), 'NOM', iocc=iocc, nbval=nbval, vect=group_detr, &
                            nbret=iret)
!              ON RECUPERE LES NUMEROS DES GROUPES A DETRUIRE
                do i = 1, nbval
                    call jenonu(jexnom(ma//group(ig), group_detr(i)), numgm)
                    if (numgm .ne. 0) then
                        zi(jgm+numgm-1) = numgm
                    end if
                end do
                AS_DEALLOCATE(vk24=group_detr)
            end do
!           ON COMPTE LE NOMBRE DE GROUPES A DETRUIRE
            nbgmde = 0
            do i = 1, nbtgp
                if (zi(jgm+i-1) .ne. 0) then
                    nbgmde = nbgmde+1
                end if
            end do
!           REACTUALISATION DE L'OBJET .GROUPEMA (OU .GROUPENO)
            nbgp = nbtgp-nbgmde
            if (nbgp .eq. 0) then
                call jedetr(ma//group(ig))
                goto 100
            end if
            call jecreo(gpptnm, 'V N K24')
            call jeecra(gpptnm, 'NOMMAX', nbgp)
            call jecrec(grp, 'V V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE', &
                        nbgp)
            do i = 1, nbtgp
                if (zi(jgm+i-1) .eq. 0) then
                    call jenuno(jexnum(ma//group(ig), i), nomgp)
                    call jecroc(jexnom(grp, nomgp))
                    call jelira(jexnom(ma//group(ig), nomgp), 'LONMAX', nbmagp)
                    call jeecra(jexnom(grp, nomgp), 'LONMAX', max(1, nbmagp))
                    call jeecra(jexnom(grp, nomgp), 'LONUTI', nbmagp)
                    call jeveuo(jexnom(grp, nomgp), 'E', j2)
                    call jeveuo(jexnom(ma//group(ig), nomgp), 'L', j1)
                    do j = 1, nbmagp
                        zi(j2+j-1) = zi(j1+j-1)
                    end do
                end if
            end do
            call jedetr(ma//group(ig))
            call jedetr(ma//ptrn(ig))
            call cpclma('&&DETGNM', ma, group(ig) (2:9), 'G')
        end if
        call jedetr('&&DETGNM.GROUP')
        call jedetr(grp)
        call jedetr(gpptnm)
100     continue
    end do
!
    call jedema()
!
end subroutine
