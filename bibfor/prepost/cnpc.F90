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

subroutine cnpc(main, macou, macsu, conneo)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utlisi.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: macou, macsu
    character(len=8) :: main
    character(len=24) :: conneo
!
!
! ----------------------------------------------------------------------
!        CONNECTIVITE DE L'ORDRE DES NOEUDS ENTRE MAILLE DE PEAU
!        ET MAILLE VOLUMIQUE CORRESPONDANTE
! ----------------------------------------------------------------------
! IN        MAIN    MAILLAGE
! IN        MACOU   NUMERO DE LA MAILLE COURANTE
! IN        MACSU   NUMERO DE LA MAILLE VOLUMIQUE CORESPONDANTE
! OUT       CONNEO  CONNECTIVIT2 DES ORDRE DES NOEUDS
! ----------------------------------------------------------------------
!
    integer(kind=8) :: inc1, inc2, nbno1, nbno2, ntrou
    integer(kind=8) :: jmacsu, jmacou, jconneo, varaux(1)

! ----------------------------------------------------------------------
!
    call jemarq()
!
    call jeveuo(jexnum(main//'.CONNEX', macsu), 'L', jmacsu)
    call jeveuo(jexnum(main//'.CONNEX', macou), 'L', jmacou)
    call jelira(jexnum(main//'.CONNEX', macou), 'LONMAX', nbno1)
    call jelira(jexnum(main//'.CONNEX', macsu), 'LONMAX', nbno2)
    call wkvect(conneo, 'V V I', nbno2, jconneo)
    do inc1 = 1, nbno2
        zi(jconneo+inc1-1) = 0
        do inc2 = 1, nbno1
            call utlisi('INTER', zi(jmacsu+inc1-1), 1, zi(jmacou+inc2-1), 1, &
                        varaux, 1, ntrou)
            if (ntrou .eq. 1) then
                zi(jconneo+inc1-1) = inc2
            end if
        end do
    end do
    call jedema()
!
! rajouter un test d'erreur
!
end subroutine
