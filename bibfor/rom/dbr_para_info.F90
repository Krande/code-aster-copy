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
subroutine dbr_para_info(cmdPara)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dbrParaInfoGreedy.h"
#include "asterfort/dbrParaInfoOrtho.h"
#include "asterfort/dbrParaInfoPod.h"
#include "asterfort/dbrParaInfoTrunc.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_ParaDBR), intent(in) :: cmdPara
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE
!
! Print informations about parameters
!
! --------------------------------------------------------------------------------------------------
!
! In  cmdPara          : datastructure for parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=16) :: operation
    aster_logical :: lReuse
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
!
! - Get parameters
!
    operation = cmdPara%operation
    lReuse = cmdPara%lReuse
!
! - Print - General
!
    if (niv .ge. 2) then
        call utmess('I', 'ROM19_6', sk=operation)
        if (lReuse) then
            call utmess('I', 'ROM19_7')
        else
            call utmess('I', 'ROM19_8')
        end if
    end if
!
! - Print / method
!
    if (operation(1:3) .eq. 'POD') then
        call dbrParaInfoPod(operation, cmdPara%paraPod)

    elseif (operation .eq. 'GLOUTON') then
        call dbrParaInfoGreedy(cmdPara%paraGreedy)

    elseif (operation .eq. 'TRONCATURE') then
        call dbrParaInfoTrunc()

    elseif (operation .eq. 'ORTHO') then
        call dbrParaInfoOrtho(cmdPara%paraOrtho)

    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
