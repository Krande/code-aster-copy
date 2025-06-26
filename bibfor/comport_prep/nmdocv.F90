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
subroutine nmdocv(keywordfact, iocc, algo_inte, keyword, l_mfront_proto, l_kit_thm, vali, valr)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: keywordfact
    integer(kind=8), intent(in) :: iocc
    character(len=16), intent(in) :: algo_inte
    character(len=14), intent(in) :: keyword
    aster_logical, intent(in) :: l_mfront_proto
    aster_logical, intent(in) :: l_kit_thm
    integer(kind=8), pointer, optional :: vali
    real(kind=8), pointer, optional :: valr
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Get and check special keywords for convergence criterion
!
! --------------------------------------------------------------------------------------------------
!
! In  keywordfact     : factor keyword to read (COMPORTEMENT)
! In  iocc            : factor keyword index in COMPORTEMENT
! In  algo_inte       : integration algorithm
! In  keyword         : keyword
! In  l_mfront_proto  : flag for a Mfront law in proto mode
! In  l_kit_thm       : flag for a law within THM kit
! Inout  vali         : pointer to the value of ITER_INTE_MAXI
! Inout  valr         : pointer to the value of RESI_INTE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: value_i, iret_i, iret_r
    real(kind=8) :: value_r
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(keyword(1:9) .eq. 'RESI_INTE' .or. keyword .eq. 'ITER_INTE_MAXI')
!
! - Get Values
!
    iret_r = 0
    iret_i = 0
    value_i = 0
    value_r = 0.d0
    if (keyword .eq. 'RESI_INTE') then
        call getvr8(keywordfact, keyword, iocc=iocc, scal=value_r, nbret=iret_r)
        ASSERT(l_mfront_proto .and. .not. l_kit_thm .or. iret_r .gt. 0)
    else if (keyword .eq. 'ITER_INTE_MAXI') then
        call getvis(keywordfact, keyword, iocc=iocc, scal=value_i, nbret=iret_i)
    end if
!
! - ITER_INTE_MAXI makes no sense for ANALYTIQUE
!
    if (algo_inte .eq. 'ANALYTIQUE') then
        if (keyword .eq. 'ITER_INTE_MAXI') then
            value_i = -value_i
        end if
    end if
!
! - Associating pointer
!
    if (iret_r .gt. 0) then
        if (.not. associated(valr)) allocate (valr)
        valr = value_r
    end if
    if (iret_i .gt. 0) then
        if (.not. associated(vali)) allocate (vali)
        vali = value_i
    end if
!
! - Checking
!
    if (keyword .eq. 'RESI_INTE') then
        if (value_r .gt. 1.0001d-6) then
            call utmess('A', 'COMPOR4_62')
        end if
    end if
!
end subroutine
