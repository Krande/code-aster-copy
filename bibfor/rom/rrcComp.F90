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
subroutine rrcComp(cmdPara)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/romFieldBuildComp.h"
#include "asterfort/rrcResultCopyParameters.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_ParaRRC), intent(in) :: cmdPara
!
! --------------------------------------------------------------------------------------------------
!
! REST_REDUIT_COMPLET - Compute
!
! Compute
!
! --------------------------------------------------------------------------------------------------
!
! In  cmdPara          : datastructure for parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iFieldBuild, nbFieldBuild, nbStore
    character(len=8) :: resultDomName, resultRomName
    type(ROM_DS_FieldBuild) :: fieldBuild
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM16_3')
    end if
!
! - Get parameters
!
    nbFieldBuild = cmdPara%nbFieldBuild
    resultDomName = cmdPara%resultDom%resultName
    resultRomName = cmdPara%resultRom%resultName
    nbStore = cmdPara%resultDom%nbStore
!
! - Compute for all fields
!
    do iFieldBuild = 1, nbFieldBuild
! ----- Current field
        fieldBuild = cmdPara%fieldBuild(iFieldBuild)

! ----- Computation for all storing index
        call romFieldBuildComp(resultDomName, resultRomName, &
                               nbStore, fieldBuild)
    end do
!
! - Copy parameters from ROM results datastructure to DOM results datastructure
!
    call rrcResultCopyParameters(cmdPara)
!
end subroutine
