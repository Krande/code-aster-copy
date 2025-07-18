! --------------------------------------------------------------------
! Copyright (C) 2016 Stefan H. Reiterer               WWW.CODE-ASTER.ORG
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
subroutine nmresx(sddisc, sderro, iter_newt)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/nmcrel.h"
#include "asterfort/nmlere.h"
#include "asterfort/utdidt.h"
#include "asterfort/getFailEvent.h"
!
    character(len=19), intent(in) :: sddisc
    character(len=24), intent(in) :: sderro
    integer(kind=8), intent(in) :: iter_newt
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Events
!
! Check if RESI_GLOB_MAXI is not too large
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for discretization
! In  sderro           : datastructure for error management (events)
! In  iter_newt        : index of current Newton iteration
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: r(1), vale_resi
    aster_logical :: l_resi_maxi
    integer(kind=8) :: nb_fail, i_fail, i_fail_acti
    integer(kind=8) :: eventType
!
! --------------------------------------------------------------------------------------------------
!
    l_resi_maxi = .false.
!
! - Index of RESI_MAXI index
!
    call utdidt('L', sddisc, 'LIST', 'NECHEC', vali_=nb_fail)
    i_fail_acti = 0
    do i_fail = 1, nb_fail
        call getFailEvent(sddisc, i_fail, eventType)
        if (eventType .eq. FAIL_EVT_RESI_MAXI) then
            i_fail_acti = i_fail
        end if
    end do
!
! - Get RESI_GLOB_MAXI
!
    call nmlere(sddisc, 'L', 'VMAXI', iter_newt, r(1))
!
! - Evaluate event
!
    if (i_fail_acti .gt. 0) then
        call utdidt('L', sddisc, 'ECHE', 'RESI_GLOB_MAXI', index_=i_fail_acti, valr_=vale_resi)
        if (r(1) .gt. vale_resi) then
            l_resi_maxi = .true.
        end if
    end if
!
! - Save event
!
    call nmcrel(sderro, 'RESI_MAXI', l_resi_maxi)
!
end subroutine
