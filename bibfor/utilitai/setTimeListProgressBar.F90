! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine setTimeListProgressBar(sddisc, nume_inst, final_)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/diinst.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
#include "asterfort/getTimeListBounds.h"
#include "jeveux.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
!
character(len=19), intent(in) :: sddisc
integer, intent(in) :: nume_inst
aster_logical, optional, intent(in) :: final_
!
! --------------------------------------------------------------------------------------------------
!
! Display progress bar for time step list
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current inst step
! In  final            : flag for final time step
!
! --------------------------------------------------------------------------------------------------
!
    integer :: perc
    real(kind=8) :: time_curr, t_ini, t_end, t_arch
    character(len=24) :: sdarch_ainf
    integer, pointer :: v_sdarch_ainf(:) => null()
    character(len=19) :: sdarch
    character(len=24) :: sdlist
    integer :: jlist, iret
!
! --------------------------------------------------------------------------------------------------
!
    time_curr = diinst(sddisc, nume_inst)
    sdarch      = sddisc(1:14)//'.ARCH'
    sdarch_ainf = sdarch(1:19)//'.AINF'
    sdlist      = sdarch(1:19)//'.LIST'
    call jeveuo(sdarch_ainf, 'L', vi = v_sdarch_ainf)
    call jeexin(sdlist, iret)


    call getTimeListBounds(sddisc, t_ini, t_end)
    if (present(final_)) then
        ASSERT(final_)
        perc = 100
        t_arch = t_end
    else
        perc = int(100.d0*(time_curr-t_ini)/(t_end-t_ini))
        if (iret.gt.0) then
            ! s'il existe une liste d'archivage
            call jeveuo(sdlist, 'L', jlist)
            t_arch = zr(jlist + v_sdarch_ainf(1)-1)
        else
            t_arch = time_curr
        endif        
    endif
    call utmess('I', 'PROGRESS_1', ni=2, vali=[perc, nume_inst],&
                                   nr=2, valr=[time_curr, t_arch])
!
end subroutine
