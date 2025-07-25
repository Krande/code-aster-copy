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

subroutine nmarpr(result, sddisc, lreuse, numder, insder, &
                  numarc)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmttch.h"
#include "asterfort/utmess.h"
    real(kind=8) :: insder
    aster_logical :: lreuse
    integer(kind=8) :: numder, numarc
    character(len=19) :: sddisc
    character(len=8) :: result
!
! ----------------------------------------------------------------------
!
! ROUTINE *_NON_LINE (ARCHIVAGE)
!
! PREMIER NUMERO A ARCHIVER
!
! ----------------------------------------------------------------------
!
!
! IN  RESULT : NOM DE LA SD RESULTAT
! IN  SDDISC : SD DISCRETISATION
! IN  NUMDER : DERNIER NUMERO ARCHIVE DANS RESULT
!               (OU 0 SI NON REENTRANT)
! IN  INSDER : DERNIER INSTANT ARCHIVE DANS RESULT
!               (R8VIDE SI NON REENTRANT)
! IN  LREUSE : .TRUE. SI CONCEPT REENTRANT
! OUT NUMARC : NUMERO DU PREMIER PAS A ARCHIVER
!
!
!
!
    character(len=24) :: tpsdit
    integer(kind=8) :: jtemps
    real(kind=8) :: valr(2), inst2
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- NOM SD_DISC
!
    tpsdit = sddisc(1:19)//'.DITR'
    call jeveuo(tpsdit, 'L', jtemps)
!
    if (lreuse) then
!
! ----- FUTUR INSTANT A ARCHIVER
!
        inst2 = zr(jtemps+2-1)
!
! ----- L'INSTANT INITIAL EST-IL SUPERIEUR AU DERNIER INSTANT ?
!
        if (inst2 .le. insder) then
            valr(1) = insder
            valr(2) = inst2
            call utmess('I', 'ARCHIVAGE_1', nr=2, valr=valr)
            call nmttch(result, inst2, numder)
            numarc = numder
        else
            numarc = numder+1
        end if
!
    else
        ASSERT(numder .eq. 0)
        numarc = 0
    end if
!
    call jedema()
!
end subroutine
