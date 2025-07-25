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
subroutine ntcrlm(listr8_sdaster, sddisc, list_inst_work)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterc/r8maem.h"
#include "asterfort/jedup1.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utdidt.h"
#include "asterfort/wkvect.h"
!
    character(len=19), intent(in) :: list_inst_work
    character(len=19), intent(in) :: sddisc
    character(len=19), intent(in) :: listr8_sdaster
!
! --------------------------------------------------------------------------------------------------
!
! THER_NON_LINE - Time discretization datastructure
!
! Create list of times and information vector from LISTR8_SDASTER
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  listr8_sdaster   : list of reals (listr8_sdaster)
! In  list_inst_work   : name of working list of time
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_inst, i_inst
    real(kind=8) :: dtmin, deltat
    character(len=8) :: list_method
    character(len=24) :: sddisc_linf
    real(kind=8), pointer :: v_sddisc_linf(:) => null()
    real(kind=8), pointer :: v_vale(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sddisc_linf = sddisc(1:19)//'.LINF'
    dtmin = r8maem()
!
! - Access to list of times
!
    call jeveuo(listr8_sdaster//'.VALE', 'L', vr=v_vale)
    call jelira(listr8_sdaster//'.VALE', 'LONMAX', nb_inst)
!
! - Minimum time between two steps
!
    do i_inst = 1, nb_inst-1
        deltat = v_vale(1+i_inst)-v_vale(i_inst)
        dtmin = min(deltat, dtmin)
    end do
!
! - Copy listr8sdaster in list of times
!
    call jedup1(listr8_sdaster(1:19)//'.VALE', 'V', list_inst_work)
!
! - Create information vector
!
    call wkvect(sddisc_linf, 'V V R', SIZE_LLINR, vr=v_sddisc_linf)
!
! - Update information vector
!
    list_method = 'MANUEL'
    call utdidt('E', sddisc, 'LIST', 'METHODE', valk_=list_method)
    call utdidt('E', sddisc, 'LIST', 'DTMIN', valr_=dtmin)
    call utdidt('E', sddisc, 'LIST', 'NBINST', vali_=nb_inst)
!
end subroutine
