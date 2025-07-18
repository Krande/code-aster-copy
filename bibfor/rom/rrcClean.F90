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
subroutine rrcClean(cmdPara)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/romFieldBuildClean.h"
#include "asterfort/romResultClean.h"
#include "asterfort/romTableClean.h"
!
    type(ROM_DS_ParaRRC), intent(inout) :: cmdPara
!
! --------------------------------------------------------------------------------------------------
!
! REST_REDUIT_COMPLET - Compute
!
! Clean datastructures
!
! --------------------------------------------------------------------------------------------------
!
! IO  cmdPara          : datastructure for parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iFieldBuild
!
! --------------------------------------------------------------------------------------------------
!
    call romResultClean()
    call romResultClean()
    do iFieldBuild = 1, cmdPara%nbFieldBuild
        call romFieldBuildClean(cmdPara%fieldBuild(iFieldBuild))
    end do
    deallocate (cmdPara%fieldBuild)
    deallocate (cmdPara%fieldName)
    call romTableClean(cmdPara%tablReduCoor)
!
end subroutine
