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
subroutine nmcrls(sddisc, listInstJv, numeInit, numeEnd, &
                  nbInstNew, dtmin)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utdidt.h"
#include "asterfort/wkvect.h"
!
    character(len=19), intent(in) :: sddisc, listInstJv
    integer(kind=8), intent(in) :: numeInit, numeEnd
    integer(kind=8), intent(out) :: nbInstNew
    real(kind=8), intent(out) :: dtmin
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Time discretization datastructure
!
! Resize list of times
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  listInstJv       : list of times from INCREMENT/LIST_INST
! In  numeInit         : index of initial time
! In  numeEnd          : index of final time
! Out nbInst           : number of time steps in list after resize
! Out dtmin            : minimum time between two steps after resize
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: pos, iInst, nbInst
    real(kind=8) :: deltat
    real(kind=8), pointer :: listInst(:) => null()
    character(len=24) :: sddiscDITRJv
    real(kind=8), pointer :: sddiscDITR(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call utdidt('L', sddisc, 'LIST', 'NBINST', vali_=nbInst)

! - Final number of time steps
    nbInstNew = (numeEnd-numeInit)+1
    ASSERT(nbInstNew .le. nbInst)

! - Acces to list of times
    call jeveuo(listInstJv, 'L', vr=listInst)

! - Create new list of time
    sddiscDITRJv = sddisc(1:19)//'.DITR'
    call wkvect(sddiscDITRJv, 'V V R', nbInstNew, vr=sddiscDITR)

! - Update new list of time
    pos = 1
    do iInst = numeInit, numeEnd
        sddiscDITR(pos) = listInst(iInst+1)
        pos = pos+1
    end do
    ASSERT(pos-1 .eq. nbInstNew)

! - New minimum time between two steps
    dtmin = r8maem()
    do iInst = 1, nbInstNew-1
        deltat = sddiscDITR(iInst+1)-sddiscDITR(iInst)
        dtmin = min(deltat, dtmin)
    end do
!
end subroutine
