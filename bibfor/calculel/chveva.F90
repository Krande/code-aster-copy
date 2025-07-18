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
subroutine chveva(nbma, ligr1, ligr2, iret)
!
! person_in_charge: jacques.pellet at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/liglma.h"
#include "asterfort/wkvect.h"
    character(len=19) :: ligr1, ligr2
    integer(kind=8) :: nbma, iret
!
! ----------------------------------------------------------------------
!
! VERIFIER LA COHERENCE DES LIGREL ENTRE DES CHAMPS PORTANT DES
! COMPOSANTES DYNAMIQUES DE TYPE VARI_ELGA
!
! ----------------------------------------------------------------------
!
!
! IN  NBMA   : NOMBRE DE MAILLES TOTALES DU MAILLAGE
! IN  LIGR1  : PREMIER LIGREL
! IN  LIGR2  : SECOND LIGREL
! OUT IRET   : 0  - SI LE SECOND LIGREL EST INCLU DANS LE PREMIER
!              -1 - SINON
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbma1, nbma2, numa, ima
    integer(kind=8) :: ipres1, ipres2
    character(len=24) :: linum1, linum2
    integer(kind=8) :: jligr1, jligr2
    character(len=24) :: linut1, linut2
    character(len=24) :: tbtrav
    integer(kind=8) :: jtrav
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- OBJETS POUR TOUT LE MAILLAGE
!
    iret = 0
    tbtrav = '&&CHVEVA.TABLE'
    call wkvect(tbtrav, 'V V I', 2*nbma, jtrav)
!
! --- EXTRACTION DU PREMIER LIGREL DE LA LISTE DES NUMEROS DE MAILLES
!
    linum1 = '&&CHVEVA.LINUM1'
    linut1 = '&&CHVEVA.LINUT1'
    call liglma(ligr1, nbma1, linum1, linut1)
    call jeveuo(linum1, 'L', jligr1)
    do ima = 1, nbma1
        numa = zi(jligr1-1+ima)
        zi(jtrav+2*(numa-1)-1+1) = 1
    end do
    call jedetr(linum1)
    call jedetr(linut1)
!
! --- EXTRACTION DU SECOND LIGREL DE LA LISTE DES NUMEROS DE MAILLES
!
    linum2 = '&&CHVEVA.LINUM2'
    linut2 = '&&CHVEVA.LINUT2'
    call liglma(ligr2, nbma2, linum2, linut2)
    call jeveuo(linum2, 'L', jligr2)
    do ima = 1, nbma2
        numa = zi(jligr2-1+ima)
        zi(jtrav+2*(numa-1)-1+2) = 1
    end do
    call jedetr(linum2)
    call jedetr(linut2)
!
! --- VERIFICATION
!
    do ima = 1, nbma
        ipres1 = zi(jtrav+2*(numa-1)-1+1)
        ipres2 = zi(jtrav+2*(numa-1)-1+2)
        if ((ipres1 .eq. 0) .and. (ipres2 .ne. 0)) then
            iret = -1
            goto 99
        end if
    end do
!
99  continue
!
    call jedetr(tbtrav)
    call jedema()
end subroutine
