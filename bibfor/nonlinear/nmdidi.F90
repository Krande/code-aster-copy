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
subroutine nmdidi(ds_inout, valinc)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/rsexch.h"
#include "asterfort/nmchex.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
!
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19), intent(in) :: valinc(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Initialization
!
! Get displacement field for DIDI loads
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_inout         : datastructure for input/output management
! In  valinc           : hat variable for algorithm fields
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, didi_nume, codret
    character(len=19) :: disp_didi, disp_prev, fieldFromResult
!
! --------------------------------------------------------------------------------------------------
!
    call nmchex(valinc, 'VALINC', 'DEPMOI', disp_prev)
    disp_didi = disp_prev
!
! - Get displacement field
!
    disp_didi = '&&CNCHAR.DIDI'
    didi_nume = ds_inout%didi_nume
    if ((didi_nume .ge. 0) .and. (ds_inout%l_stin_evol)) then
        call rsexch(' ', ds_inout%stin_evol, 'DEPL', didi_nume, fieldFromResult, iret)
        if (iret .ne. 0) then
            call utmess('F', 'MECANONLINE5_20', sk=ds_inout%stin_evol)
        end if
    else
        fieldFromResult = disp_prev
    end if
!
    call vtcopy(fieldFromResult, disp_didi, codret)
    if (codret .ne. 0) then
        call utmess("F", "FIELD0_8")
    end if
!
end subroutine
