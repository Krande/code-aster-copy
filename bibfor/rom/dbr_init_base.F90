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
subroutine dbr_init_base(cmdPara)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dbrInitBaseGreedy.h"
#include "asterfort/dbrInitBaseOrtho.h"
#include "asterfort/dbrInitBasePod.h"
#include "asterfort/dbrInitBaseTrunc.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_ParaDBR), intent(inout) :: cmdPara
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE
!
! Initializations for base
!
! --------------------------------------------------------------------------------------------------
!
! IO  cmdPara          : datastructure for parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM19_4')
    end if
!
    if (cmdPara%operation(1:3) .eq. 'POD') then
        call dbrInitBasePod(cmdPara%resultOut%resultName, cmdPara%paraPod, &
                            cmdPara%lReuse, cmdPara%base)
    elseif (cmdPara%operation .eq. 'GLOUTON') then
        call dbrInitBaseGreedy(cmdPara%resultOut%resultName, cmdPara%paraGreedy, &
                               cmdPara%base)
    elseif (cmdPara%operation .eq. 'TRONCATURE') then
        call dbrInitBaseTrunc(cmdPara%resultOut%resultName, cmdPara%paraTrunc, &
                              cmdPara%lReuse, cmdPara%base)
    elseif (cmdPara%operation .eq. 'ORTHO') then
        call dbrInitBaseOrtho(cmdPara%resultOut%resultName, cmdPara%paraOrtho, &
                              cmdPara%lReuse, cmdPara%base)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
