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
subroutine dbrInitBaseOrtho(resultName, paraOrtho, lReuse, base)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/romBaseCreate.h"
#include "asterfort/romBaseDSCopy.h"
#include "asterfort/romBaseGetInfo.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: resultName
    type(ROM_DS_ParaDBR_Ortho), intent(inout) :: paraOrtho
    aster_logical, intent(in) :: lReuse
    type(ROM_DS_Empi), intent(inout) :: base
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE
!
! Initializations for base - For orthogonalization
!
! --------------------------------------------------------------------------------------------------
!
! In  resultName       : name of results datastructure to save base
! IO  paraOrtho        : datastructure for parameters (orthogonalization)
! In  lReuse           : .true. if reuse
! IO  base             : base
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM18_22')
    end if
!
! - Get informations about base to orthogonalize
!
    if (lReuse) then
        if (niv .ge. 2) then
            call utmess('I', 'ROM18_23')
        end if
        call romBaseGetInfo(resultName, base)
    else
        if (niv .ge. 2) then
            call utmess('I', 'ROM18_24')
        end if
        call romBaseGetInfo(paraOrtho%baseInitName, paraOrtho%baseInit)
    end if
!
! - Create base (if necessary)
!
    if (.not. lReuse) then
        if (niv .ge. 2) then
            call utmess('I', 'ROM18_25')
        end if
        call romBaseDSCopy(paraOrtho%baseInit, resultName, base)
        call romBaseCreate(base, paraOrtho%baseInit%nbMode)
    end if
!
! - If reuse: check that name is the name between output result end keyword BASE
!
    if (lReuse) then
        if (paraOrtho%baseInitName .ne. resultName) then
            call utmess('F', 'ROM18_21')
        end if
    end if
!
end subroutine
