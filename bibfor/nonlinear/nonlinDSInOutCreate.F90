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
subroutine nonlinDSInOutCreate(phenom, ds_inout)
!
    use NonLin_Datastructure_type
    use listLoad_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/CreateInOutDS_M.h"
#include "asterfort/CreateInOutDS_T.h"
!
    character(len=4), intent(in) :: phenom
    type(NL_DS_InOut), intent(out) :: ds_inout
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output management
!
! Create input/output datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  phenom           : phenomenon (MECA/THER/THNL) /VIBR
! Out ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: listLoadResu
!
! --------------------------------------------------------------------------------------------------
!
    ds_inout%result = ' '
    ds_inout%l_temp_nonl = ASTER_FALSE
    ds_inout%stin_evol = ' '
    ds_inout%l_stin_evol = ASTER_FALSE
    ds_inout%l_state_init = ASTER_FALSE
    ds_inout%l_reuse = ASTER_FALSE
    ds_inout%didi_nume = -1
    ds_inout%criterion = ' '
    ds_inout%precision = r8vide()
    ds_inout%user_time = r8vide()
    ds_inout%user_nume = 0
    ds_inout%stin_time = r8vide()
    ds_inout%l_stin_time = ASTER_FALSE
    ds_inout%l_user_time = ASTER_FALSE
    ds_inout%l_user_nume = ASTER_FALSE
    ds_inout%init_time = r8vide()
    ds_inout%init_nume = -1
    ds_inout%l_init_stat = ASTER_FALSE
    ds_inout%l_init_vale = ASTER_FALSE
    ds_inout%temp_init = r8vide()

! - Generate name of list of loads saved in results datastructure
    call nameListLoad(listLoadResu)
    ds_inout%listLoadResu = listLoadResu
!
! - Specific parameters
!
    if (phenom .eq. 'MECA') then
        call CreateInOutDS_M(ds_inout)
    elseif (phenom .eq. 'THER') then
        call CreateInOutDS_T(ds_inout, ASTER_FALSE)
    elseif (phenom .eq. 'THNL') then
        call CreateInOutDS_T(ds_inout, ASTER_TRUE)
    elseif (phenom .eq. 'VIBR') then
        call CreateInOutDS_M(ds_inout)
    else
        ASSERT(.false.)
    end if
!
end subroutine
