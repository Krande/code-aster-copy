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

subroutine acevtr(nomo, ityp, noms, itab, nn, idim)
!.======================================================================
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/testli.h"
#include "asterfort/utmess.h"
!
    character(len=24) :: noms(*)
    character(len=8) :: nomo
    integer(kind=8) :: ityp, nn, idim, itab(*)
!
!   ARGUMENT        E/S  TYPE         ROLE
!.========================= DEBUT DES DECLARATIONS ====================
! -----  VARIABLES LOCALES
    integer(kind=8) :: repi
    character(len=16) :: nomte, nomodl, chaine
    character(len=19) :: nolig
    character(len=24) :: repk
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ialiel, ierr, igrel, ima, iret, itypel
    integer(kind=8) :: kma, kmai, nbgrel, nel
!-----------------------------------------------------------------------
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    repi = 0
    if (idim .eq. 2) then
        chaine = '2D_DIS_TR'
    else
        chaine = 'DIS_TR'
    end if
!
    call dismoi('NOM_LIGREL', nomo, 'MODELE', repk=nolig)
    call jeexin(nolig//'.LIEL', iret)
    if (iret .ne. 0) then
        call jelira(nolig//'.LIEL', 'NUTIOC', nbgrel)
! le nombre de grels du LIGREL du modele est nul.
        ASSERT(nbgrel .gt. 0)
        nomodl = ' '
        ierr = 0
        do igrel = 1, nbgrel
            call jeveuo(jexnum(nolig//'.LIEL', igrel), 'L', ialiel)
            call jelira(jexnum(nolig//'.LIEL', igrel), 'LONMAX', nel)
            itypel = zi(ialiel-1+nel)
            call jenuno(jexnum('&CATA.TE.NOMTE', itypel), nomte)
            call dismoi('MODELISATION', nomte, 'TYPE_ELEM', repk=repk)
            nomodl = repk(1:16)
            if (nomodl .ne. chaine) then
                if (ityp .eq. 0) then
                    ierr = 1
                    kmai = zi(ialiel)
                    goto 20
                else
                    do kma = 1, nn
                        if (ityp .eq. 1) then
                            ima = char8_to_int(noms(kma))
                        else
                            ima = itab(kma)
                        end if
                        call testli(ima, zi(ialiel), nel-1, kmai, ierr)
                        if (ierr .eq. 1) goto 20
                    end do
                end if
            end if
        end do
    end if
20  continue
!     IF (IERR.EQ.1)  WRITE(*,*) 'KMAI',KMAI,'IGREL',IGREL,
!    .       'NOMODL',NOMODL,'CHAINE',CHAINE
    if (ierr .eq. 1) then
        call utmess('F', 'DISCRETS_9')
    end if
    call jedema()
end subroutine
