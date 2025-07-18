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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine calcGetData(table_new, table_old, &
                       nb_option, list_option, &
                       nume_inst, list_inst, &
                       phenom, l_pred)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
!
    character(len=8), intent(out) :: table_new
    character(len=8), intent(out) :: table_old
    integer(kind=8), intent(out) :: nb_option
    character(len=16), intent(out) :: list_option(:)
    integer(kind=8), intent(out) :: nume_inst
    character(len=19), intent(out) :: list_inst
    character(len=16), intent(out) :: phenom
    aster_logical, intent(out) :: l_pred
!
! --------------------------------------------------------------------------------------------------
!
! Command CALCUL
!
! Get commons data
!
! --------------------------------------------------------------------------------------------------
!
! Out table_new        : name of created table
! Out table_old        : name of old table
! Out nb_option        : number of options to compute
! Out list_option      : list of options to compute
! Out nume_inst        : index of current step time
! Out list_inst        : list of step time
! Out phenom           : phenomenon (MECANIQUE/THERMIQUE/ACOUSTIQUE)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: k24dummy
    character(len=16) :: pred
    integer(kind=8) :: nocc
!
! --------------------------------------------------------------------------------------------------
!
    table_new = ' '
    table_old = ' '
    nb_option = 0
    list_option(:) = ' '
    nume_inst = 0
    list_inst = ' '
    phenom = ' '
!
! - Name of created table
!
    call getres(table_new, k24dummy, k24dummy)
!
! - Name of reused table
!
    call getvid(' ', 'TABLE', nbval=0, nbret=nocc)
    if (nocc .eq. 0) then
        table_old = ' '
    else
        call getvid(' ', 'TABLE', nbval=1, scal=table_old)
        if (table_old .ne. table_new) then
            call utmess('F', 'CALCUL1_3')
        end if
    end if
!
! - Options
!
    call getvtx(' ', 'OPTION', nbval=6, vect=list_option, nbret=nb_option)

    call getvtx(' ', 'PHASE', scal=pred, nbret=nocc)

    if (nocc .eq. 1 .and. pred .eq. 'PREDICTION') then
        l_pred = ASTER_TRUE
    else
        l_pred = ASTER_FALSE
    end if
!
! - Phenomen
!
    call getvtx(' ', 'PHENOMENE', scal=phenom, nbret=nocc)
    ASSERT(nocc .eq. 1)
!
! - Get current time
!
    call getvis('INCREMENT', 'NUME_ORDRE', iocc=1, scal=nume_inst)
    call getvid('INCREMENT', 'LIST_INST', iocc=1, scal=list_inst)
!
end subroutine
