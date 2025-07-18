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
subroutine romBaseGetInfo(resultName, base)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/romModeParaRead.h"
#include "asterfort/rs_getfirst.h"
#include "asterfort/rs_get_liststore.h"
#include "asterfort/rsexch.h"
#include "asterfort/romFieldGetInfo.h"
!
    character(len=8), intent(in)     :: resultName
    type(ROM_DS_Empi), intent(inout) :: base
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Base management
!
! Get informations about base
!
! --------------------------------------------------------------------------------------------------
!
! In  resultName       : name of results datastructures for base
! IO  base             : base
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: numeModeRefe = 1
    integer(kind=8) :: iret, numeModeFirst, numeSlice, nbSnap, nbMode
    character(len=8)  :: model, lineicAxis, baseType
    character(len=24) :: lineicSect, fieldRefe, modeSymbName
    type(ROM_DS_Field) :: mode
!
! --------------------------------------------------------------------------------------------------
!
    nbMode = 0
    model = ' '
    lineicAxis = ' '
    baseType = ' '
    lineicSect = ' '
    fieldRefe = ' '
    modeSymbName = ' '
!
! - Number of modes
!
    call rs_get_liststore(resultName, nbMode)
!
! - Get main parameters in empiric result
!
    call romModeParaRead(resultName, numeModeRefe, &
                         model_=model, &
                         modeSymbName_=modeSymbName, &
                         numeSlice_=numeSlice, &
                         nbSnap_=nbSnap)
    if (numeSlice .eq. 0) then
        baseType = '3D'
    else
        baseType = 'LINEIQUE'
    end if
!
! - Get _representative_ field in empiric result
!
    call rs_getfirst(resultName, numeModeFirst)
    call rsexch(' ', resultName, modeSymbName, numeModeFirst, fieldRefe, iret)
    ASSERT(iret .eq. 0)
!
! - Get informations from mode
!
    call romFieldGetInfo(model, modeSymbName, fieldRefe, mode)
!
! - Save informations about empiric modes
!
    base%resultName = resultName
    base%baseType = baseType
    base%lineicAxis = lineicAxis
    base%lineicSect = lineicSect
    base%nbMode = nbMode
    base%nbSnap = nbSnap
    base%mode = mode
!
end subroutine
