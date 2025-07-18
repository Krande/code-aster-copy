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

subroutine chpve2(nomch, nbtyp, tabtyp, ier)
!
!     VERIFICATIONS DE LA GRANDEUR ET DE LA LOCALISATION DES CHAMPS.
!
!  IN  NOCHAM : NOM DU CHAMP
!  IN  NBTYP  : DIMENSION DE TABTYP
!  IN  TABTYP : TABLEAU CONTENANT LES TYPES DE CHAMPS ACCEPTABLES.
!               UN ELEMENT DE TABTYP EST DE LA FORME : LOC#GD
!               OU : LOC = ELNO/ELGA/ELEM/ELXX/CART
!                    GD  = GRANDEUR
!  OUT   IERD  : CODE RETOUR  (0--> OK, 1--> PB )
! ======================================================================
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ier, nbtyp
    character(len=*) :: tabtyp(nbtyp), nomch
!
    integer(kind=8) :: lc, i, j
    character(len=19) :: noch
    character(len=4) :: loch, tych
    character(len=8) :: gdch, nomgd, blan8
    character(len=11) :: chaine
    character(len=24) :: valk
!
    call jemarq()
!
    ier = 1
    noch = nomch
    blan8 = '        '
    nomgd = blan8
    gdch = blan8
    do i = 1, nbtyp
        lc = len(tabtyp(i))
        ASSERT(lc .ge. 11)
        chaine = tabtyp(i) (1:11)
!
        do j = 1, lc
            if (chaine(j:j) .eq. '#') then
                loch = chaine(1:j)
                gdch = chaine(j+1:11)
                goto 30
            end if
        end do
30      continue
!
        call dismoi('TYPE_CHAMP', noch, 'CHAMP', repk=tych)
        call dismoi('NOM_GD', noch, 'CHAMP', repk=nomgd)
!
        if ((loch(3:4) .ne. 'XX' .and. loch .eq. tych) .or. &
            (loch(3:4) .eq. 'XX' .and. loch(1:2) .eq. tych(1:2))) then
            if (gdch(1:6) .eq. nomgd(1:6)) then
                ier = 0
                goto 40
            end if
        end if
    end do
40  continue
!
    if (ier .ne. 0) then
        valk = tych//'_'//nomgd
        call utmess('F', 'UTILITAI5_97', sk=valk)
    end if
!
    call jedema()
!
end subroutine
