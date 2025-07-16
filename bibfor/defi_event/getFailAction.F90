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
subroutine getFailAction(sddisc, failIndx, actionType)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/jeveuo.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: failIndx
    integer(kind=8), intent(out) :: actionType
!
! --------------------------------------------------------------------------------------------------
!
! Event management
!
! Get action type for failure
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : name of datastructure for time discretization
! In  failIndx         : current index for ECHEC keyword
! Out actionType       : type of action
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sddiscEEVRJv
    real(kind=8), pointer :: sddiscEEVR(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sddiscEEVRJv = sddisc(1:19)//'.EEVR'
    call jeveuo(sddiscEEVRJv, 'L', vr=sddiscEEVR)
    actionType = nint(sddiscEEVR(SIZE_LEEVR*(failIndx-1)+2))
!
end subroutine
