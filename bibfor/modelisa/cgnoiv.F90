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

subroutine cgnoiv(iocc, nomaz, lisnoz, nbno)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: iocc, nbno
    character(len=*) :: nomaz, lisnoz
! person_in_charge: jacques.pellet at edf.fr
!
!       CGNOIV -- TRAITEMENT DE L'OPTION INTERVALLE_VALE
!                 DU MOT FACTEUR CREA_GROUP_NO DE
!                 LA COMMANDE DEFI_GROUP
!
! -------------------------------------------------------
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  LISNOZ        - JXVAR - K24  - : NOM DE LA LISTE DE NOEUDS RETENUS
!  NBNO          - OUT   -  I   - : LONGUEUR DE CETTE LISTE
! -------------------------------------------------------
!
    integer(kind=8) :: n1, k, nbnot, ino, ncmp
    integer(kind=8) :: jlisno, jcn2v, jcn2l
    character(len=3) :: tsca
    character(len=8) :: nocmp, noma, ma1, nomgd
    character(len=16) :: motfac
    character(len=19) :: cham19, cns1, cns2
    character(len=24) :: lisnoi, valk(5)
    real(kind=8) :: valr(2), vmin, vmax, v1
    character(len=8), pointer :: cnsk(:) => null()
    integer(kind=8), pointer :: lisno(:) => null()
!     -----------------------------------------------------------------
!
    call jemarq()
!
    motfac = 'CREA_GROUP_NO'
    lisnoi = lisnoz
    noma = nomaz
    nbno = 0
!
    call getvid(motfac, 'CHAM_GD', iocc=iocc, scal=cham19, nbret=n1)
    call getvtx(motfac, 'NOM_CMP', iocc=iocc, scal=nocmp, nbret=n1)
    call getvr8(motfac, 'VALE', iocc=iocc, nbval=2, vect=valr(1), &
                nbret=n1)
    ASSERT(n1 .eq. 2)
    vmin = valr(1)
    vmax = valr(2)
    if (vmin .gt. vmax) goto 30
!
    call dismoi('NOM_MAILLA', cham19, 'CHAMP', repk=ma1)
    if (noma .ne. ma1) then
        valk(1) = cham19
        valk(2) = noma
        call utmess('F', 'CALCULEL2_50', nk=2, valk=valk)
    end if
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbnot)
!
    cns1 = '&&CGNOIV.CNS1'
    cns2 = '&&CGNOIV.CNS2'
    call cnocns(cham19, 'V', cns1)
    ncmp = 1
    call cnsred(cns1, 0, [0], ncmp, nocmp, &
                'V', cns2)
!
    call jeveuo(cns2//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(cns2//'.CNSV', 'L', jcn2v)
    call jeveuo(cns2//'.CNSL', 'L', jcn2l)
    nomgd = cnsk(2)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    ASSERT(tsca .eq. 'R' .or. tsca .eq. 'I')
!
    AS_ALLOCATE(vi=lisno, size=nbnot)
    do ino = 1, nbnot
        if (.not. zl(jcn2l-1+(ino-1)*ncmp+1)) goto 10
        if (tsca .eq. 'R') then
            v1 = zr(jcn2v-1+(ino-1)*ncmp+1)
        else if (tsca .eq. 'I') then
            v1 = zi(jcn2v-1+(ino-1)*ncmp+1)
        end if
        if (v1 .ge. vmin .and. v1 .le. vmax) then
            nbno = nbno+1
            lisno(nbno) = ino
        end if
10      continue
    end do
!
!
!
! --- ALLOCATION DU VECTEUR DES NUMEROS DES MAILLES RETENUES
!     --------------------------------------------------------
    call wkvect(lisnoi, 'V V I', max(nbno, 1), jlisno)
    do k = 1, nbno
        zi(jlisno-1+k) = lisno(k)
    end do
!
!
! --- MENAGE :
    AS_DEALLOCATE(vi=lisno)
    call detrsd('CHAM_NO_S', cns1)
    call detrsd('CHAM_NO_S', cns2)
!
30  continue
    call jedema()
!
end subroutine
