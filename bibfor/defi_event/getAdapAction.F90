! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine getAdapAction(sddisc, iAdap, actionType)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jeveuo.h"
#include "event_def.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: iAdap
    integer(kind=8), intent(out) :: actionType
!
! --------------------------------------------------------------------------------------------------
!
! Event management
!
! Get action type for time adaptation
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : name of datastructure for time discretization
! In  iAdap            : current index for ADAPTATION keyword
! Out actionType       : type of action
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: addiscAptrJv
    real(kind=8), pointer :: sddiscAptr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    addiscAptrJv = sddisc(1:19)//'.ATPR'
    call jeveuo(addiscAptrJv, 'L', vr=sddiscAptr)
    actionType = nint(sddiscAptr(SIZE_LATPR*(iAdap-1)+1))
!
end subroutine
