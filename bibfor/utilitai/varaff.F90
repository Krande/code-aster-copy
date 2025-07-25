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

subroutine varaff(noma, gran, base, ceselz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    character(len=1) :: base
    character(len=8) :: noma, gran
    character(len=*) :: ceselz
! BUT :
!  - TRAITER L'OPTION 'AFFE' DE LA COMMANDE CREA_CHAMP POUR VARI_R
!  - CREER LE CHAM_ELEM_S / ELEM  (CESELZ)
!-----------------------------------------------------------------------
    integer(kind=8) :: nocc, nbtou, n1
    integer(kind=8) :: iad, k, iocc, nbmail
    character(len=8) :: kbid, typmcl(2)
    character(len=6) :: knuva
    character(len=16) :: motclf, motcls(2)
    character(len=19) :: ceselm
    character(len=24) :: mesmai
    integer(kind=8) :: nvarmx, jcesd, jcesl, jlvavx
    integer(kind=8) :: jmesma, kvari, n2, numa, nuva, nuvamx, nbmato
    parameter(nvarmx=10000)
    aster_logical :: ltou
    character(len=8), pointer :: lnova(:) => null()
    character(len=8), pointer :: lnovx(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
!   ------------------------------------------------------------------
    call jemarq()
!
    if (noma .eq. ' ') then
        call utmess('F', 'UTILITAI_10')
    end if
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbmato)
!
    ASSERT(gran .eq. 'VARI_R')
    AS_ALLOCATE(vk8=lnovx, size=nvarmx)
    call wkvect('&&VARAFF.LVAVX', 'V V R', nvarmx, jlvavx)
!
    motclf = 'AFFE'
    call getfac(motclf, nocc)
!
    mesmai = '&&VARAFF.MES_MAILLES'
    motcls(1) = 'GROUP_MA'
    motcls(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
!
!     0- CALCUL DU PLUS GRAND NUMERO DE VARI UTILISE (NUVAMX):
!     --------------------------------------------------------
    nuvamx = 0
    do iocc = 1, nocc
        call getvtx(motclf, 'NOM_CMP', iocc=iocc, nbval=nvarmx, vect=lnovx, &
                    nbret=n1)
        ASSERT(n1 .gt. 0)
        do k = 1, n1
            ASSERT(lnovx(k) (1:1) .eq. 'V')
            read (lnovx(k) (2:8), '(I7)') nuva
            nuvamx = max(nuvamx, nuva)
        end do
    end do
    ASSERT(nuvamx .gt. 0)
!
!
!     1- ALLOCATION DE CESELM
!     --------------------------------------------
    AS_ALLOCATE(vk8=lnova, size=nuvamx)
    do k = 1, nuvamx
        call codent(k, 'G', knuva)
        lnova(k) = 'V'//knuva
    end do
    ceselm = ceselz
!     -- REMARQUE : LES CMPS SERONT DANS L'ORDRE V1,V2,...
    call cescre(base, ceselm, 'ELEM', noma, 'VARI_R', &
                nuvamx, lnova, [0], [-1], [-nuvamx])
!
    call jeveuo(ceselm//'.CESD', 'L', jcesd)
    call jeveuo(ceselm//'.CESV', 'E', vr=cesv)
    call jeveuo(ceselm//'.CESL', 'E', jcesl)
!
!
!     2- BOUCLE SUR LES OCCURENCES DU MOT CLE AFFE
!     --------------------------------------------
    do iocc = 1, nocc
!
        call getvtx(motclf, 'NOEUD', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call utmess('F', 'UTILITAI_12')
        end if
        call getvtx(motclf, 'GROUP_NO', iocc=iocc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call utmess('F', 'UTILITAI_13')
        end if
!
        call getvtx(motclf, 'NOM_CMP', iocc=iocc, nbval=nvarmx, vect=lnovx, &
                    nbret=n1)
        call getvr8(motclf, 'VALE', iocc=iocc, nbval=nvarmx, vect=zr(jlvavx), &
                    nbret=n2)
        ASSERT(n1 .eq. n2)
!
!
!
        call getvtx(motclf, 'TOUT', iocc=iocc, scal=kbid, nbret=nbtou)
        if (nbtou .eq. 1) then
            ltou = .true.
            nbmail = nbmato
        else
            ltou = .false.
            call reliem(' ', noma, 'NU_MAILLE', motclf, iocc, &
                        2, motcls, typmcl, mesmai, nbmail)
            ASSERT(nbmail .gt. 0)
            call jeveuo(mesmai, 'L', jmesma)
        end if
!
        do kvari = 1, n1
            read (lnovx(kvari) (2:8), '(I7)') nuva
            ASSERT(nuva .gt. 0 .and. nuva .le. nuvamx)
            do k = 1, nbmail
                if (ltou) then
                    numa = k
                else
                    numa = zi(jmesma-1+k)
                end if
                ASSERT(numa .gt. 0 .and. numa .le. nbmato)
!
                call cesexi('C', jcesd, jcesl, numa, 1, &
                            1, nuva, iad)
                ASSERT(iad .ne. 0)
!               -- On peut vouloir surcharger :
                iad = abs(iad)
!
!               -- RECOPIE DE LA VALEUR:
                zl(jcesl-1+iad) = .true.
                cesv(iad) = zr(jlvavx-1+kvari)
            end do
        end do
!
        call jedetr(mesmai)
    end do
!
!
    AS_DEALLOCATE(vk8=lnovx)
    AS_DEALLOCATE(vk8=lnova)
    call jedema()
end subroutine
