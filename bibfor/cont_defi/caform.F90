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

subroutine caform(contForm)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cazouu.h"
#include "asterfort/getvtx.h"
#include "Contact_type.h"
!
    integer(kind=8), intent(out) :: contForm
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Get contact formulation
!
! --------------------------------------------------------------------------------------------------
!
! Out contForm        : formulation of contact
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: zoneKeyword = "ZONE"
    character(len=16) :: formulation, algoCont
    integer(kind=8) :: noc, nbContZone
!
! --------------------------------------------------------------------------------------------------
!
    contForm = CONT_FORM_UNDEF

! - Contact formulation
    call getvtx(' ', 'FORMULATION', scal=formulation, nbret=noc)
    ASSERT(noc .ne. 0)
!
    if (formulation .eq. 'DISCRETE') then
        contForm = CONT_FORM_DISC
    else if (formulation .eq. 'CONTINUE') then
        call getvtx(zoneKeyword, 'ALGO_CONT', iocc=1, scal=algoCont)
        if (algoCont .eq. 'LAC') then
            call getfac(zoneKeyword, nbContZone)
            call cazouu(zoneKeyword, nbContZone, 'ALGO_CONT', 'T')
            contForm = CONT_FORM_LAC
        else
            contForm = CONT_FORM_CONT
        end if
    else if (formulation .eq. 'LIAISON_UNIL') then
        contForm = CONT_FORM_UNIL
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
