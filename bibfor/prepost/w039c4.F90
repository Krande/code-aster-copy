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

subroutine w039c4(carte, ifi, form)
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesred.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/ircame.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
    character(len=*) :: carte, form
    integer(kind=8) :: ifi
!
! ----------------------------------------------------------------------
!
!                  IMPRIMER UNE "CARTE" AU FORMAT MED
!
!     * FORMAT MED
!     * LA CARTE DOIT EXISTER
!     * LA CARTE NE DOIT CONTENIR QUE DES REELS
!     * TOUTES LES COMPOSANTES SONT IMPRIMEES (NON AFFECTEES ==> 0)
!
! ----------------------------------------------------------------------
!
!
    integer(kind=8) :: iret, nugd, n1, jnocmp, nbCmpDyna
    integer(kind=8) :: jcesk, jcesd, jcesc, jcesv, jcesl
    character(len=8) :: typech, tsca, nomgd, k8bid
    character(len=19) :: cart1, chels1, chels2
    character(len=16) :: field_type
    character(len=64) :: nommed
    integer(kind=8), pointer :: desc(:) => null()
! ----------------------------------------------------------------------

    call jemarq()

!   --- si ce n'est pas au format med
    if (form .ne. 'MED') goto 999

!   --- si la carte n'existe pas
    call exisd('CARTE', carte, iret)
    if (iret .eq. 0) goto 999

    cart1 = carte
!   --- que des reels
    call jeveuo(cart1//'.DESC', 'L', vi=desc)
    nugd = desc(1)
    call jenuno(jexnum('&CATA.GD.NOMGD', nugd), nomgd)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    ASSERT(tsca .eq. 'R')

! --- ON TRANSFORME LA CARTE EN UN CHAM_ELEM_S
    chels1 = '&&W039C4.CHELS1'
    call carces(cart1, 'ELEM', ' ', 'V', chels1, &
                'A', iret)

!     -- LE FORMAT MED REFUSE DE TRAITER PLUS DE 80 CMPS :
    call jelira(jexnum('&CATA.GD.NOMCMP', nugd), 'LONMAX', n1)
    if (n1 .gt. 80) then
!       -- ON NE GARDE QUE LES 80 PREMIERES :
        call utmess('A', 'CALCULEL4_24', sk=cart1)
        chels2 = '&&W039C4.CHELS2'
        call copisd('CHAM_ELEM_S', 'V', chels1, chels2)
        call jeveuo(jexnum('&CATA.GD.NOMCMP', nugd), 'L', jnocmp)
        call cesred(chels2, 0, [0], 80, zk8(jnocmp), &
                    'V', chels1)
        call detrsd('CHAM_ELEM_S', chels2)
    end if
!
!
! --- POUR AVOIR UN NOM MED PROCHE DE CELUI DE LA CARTE.
!        PAS DE '_' DEJA UTILISE PAR W039C1
    nommed = carte
    nommed(9:9) = '#'
    typech = 'ELEM'
!
! --- ON RECUPERE LES OBJETS
    call jeveuo(chels1//'.CESK', 'L', jcesk)
    call jeveuo(chels1//'.CESD', 'L', jcesd)
    call jeveuo(chels1//'.CESC', 'L', jcesc)
    call jeveuo(chels1//'.CESV', 'L', jcesv)
    call jeveuo(chels1//'.CESL', 'L', jcesl)
!
! --- ECRITURE DES CHAMPS AU FORMAT MED
    k8bid = ' '
    field_type = 'Unknown'
    call ircame(ifi, nommed, chels1, typech, k8bid, &
                0, k8bid, k8bid, k8bid, 0, &
                0.0d0, 0, jcesk, jcesd, jcesc, &
                jcesv, jcesl, 0, [0], k8bid, k8bid, &
                field_type, nbCmpDyna, .false._1, iret)
!
    ASSERT(iret .eq. 0)
!
    call detrsd('CHAM_ELEM_S', chels1)
!
999 continue
    call jedema()
end subroutine
