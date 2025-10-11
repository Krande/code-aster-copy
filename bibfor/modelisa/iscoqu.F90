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
subroutine iscoqu(nomo, numail, lcoque)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
    character(len=8) :: nomo
    integer(kind=8) :: numail
    aster_logical :: lcoque
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (TOUTES METHODES - LECTURE DONNEES)
!
! DETECTE SI UN ELEMENT EST DE TYPE COQUE_3D
!
! ----------------------------------------------------------------------
!
!
! IN  NOMO   : NOM DU MODELE
! IN  NUMAIL : NUMERO ABSOLU DE LA MAILLE
! OUT LCOQUE : .TRUE. SI COQUE_3D
!
!
!
!
    integer(kind=8) :: iret, igrel, iel
    integer(kind=8) :: ialiel, itypel
    integer(kind=8) :: nbgrel, nel, numai2
    character(len=8) :: nomte
    character(len=19) :: ligrmo
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    lcoque = .false.
!
! --- LIGREL DU MODELE
!
    call dismoi('NOM_LIGREL', nomo, 'MODELE', repk=ligrmo)
    call jeexin(ligrmo//'.LIEL', iret)
    if (iret .eq. 0) then
        ASSERT(.false.)
    end if
!
! --- NOMBRE DE GREL
!
    call jelira(ligrmo(1:19)//'.LIEL', 'NUTIOC', nbgrel)
!
! --- BOUCLE SUR LES GREL
!
    do igrel = 1, nbgrel
!
! --- TYPE DU GREL COURANT
!
        call jeveuo(jexnum(ligrmo(1:19)//'.LIEL', igrel), 'L', ialiel)
        call jelira(jexnum(ligrmo(1:19)//'.LIEL', igrel), 'LONMAX', nel)
        itypel = zi(ialiel-1+nel)
        call jenuno(jexnum('&CATA.TE.NOMTE', itypel), nomte)
!
! --- CAS DES COQUES_3D
!
        if ((nomte .eq. 'MEC3QU9H') .or. (nomte .eq. 'MEC3TR7H')) then
!
! --- BOUCLE DANS LE GREL
!
            do iel = 1, nel-1
                numai2 = zi(ialiel-1+iel)
                if (numai2 .eq. numail) then
                    lcoque = .true.
                    goto 40
                end if
            end do
        end if
40      continue
    end do
!
    call jedema()
end subroutine
