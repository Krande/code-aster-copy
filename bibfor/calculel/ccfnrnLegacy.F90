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
subroutine ccfnrnLegacy(option, postComp)
!
    use postComp_type
    use postComp_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcreb.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    type(POST_COMP), intent(inout) :: postComp
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP
!
! Compute FORC_NODA and REAC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option            : option to compute
! IO  postComp          : general type for post-processing
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iStore, numeStore, iret
    character(len=24) :: fieldOut
    character(len=14) :: numeDof
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

    ASSERT(postComp%lPostNoda)

! - Init parameters to compute option
    call initPara(option, postComp)

! - Loop on storing indexes
    do iStore = 1, postComp%postCompResu%nbStore
        call jemarq()

! ----- Save the previous time and the previous generalized stresses
        postComp%postCompPara%timePrev = postComp%postCompPara%time
        postComp%postCompFields%sigmPrev = postComp%postCompFields%sigm

! ----- Current storing index
        numeStore = postComp%postCompResu%listStore(iStore)

! ----- Get parameters to compute option
        call getPara(option, numeStore, postComp)

! ----- Get input fields to compute option
        call getInputFields(option, numeStore, postComp)

! ----- Get numbering
        call getNumbering(option, postComp, numeDof)

! ----- Get access of output field
        call rsexch(' ', postComp%postCompResu%resultOut, option, numeStore, fieldOut, iret)
        call jeexin(fieldOut(1:19)//'.REFE', iret)
        if (iret .ne. 0) then
            call utmess('I', 'POSTCOMP1_11', sk=option, si=numeStore)
            call detrsd('CHAM_NO', fieldOut)
        end if

! ----- Create output field
        if (postComp%lCplx) then
            call vtcreb(fieldOut, 'G', 'C', nume_ddlz=numeDof)
        else
            call vtcreb(fieldOut, 'G', 'R', nume_ddlz=numeDof)
        end if

! ----- Compute FORC_NODA
        if (postComp%postCompNoda%lForcNoda) then
            call compForcNoda(postComp, numeDof, fieldOut)
        end if

! ----- Compute FEXT_NODA
        if (postComp%postCompNoda%lFextNoda) then
            call compFextNoda(postComp, option, numeDof, fieldOut)
        end if

! ----- Compute M_GAMMA
        if (postComp%postCompNoda%lMGamma) then
            call compMGamma(postComp, option, numeStore, numeDof, fieldOut)
        end if

! ----- Store field in result
        call rsnoch(postComp%postCompResu%resultOut, option, numeStore)
        call jedema()
    end do
!
    call jedema()
end subroutine
