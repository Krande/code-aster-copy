! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine pmdocr(carcri)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/carc_info.h"
#include "asterfort/carc_chck.h"
#include "asterfort/carc_read.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/setBehaviourParaValue.h"
#include "asterfort/Behaviour_type.h"
!
real(kind=8), intent(out) :: carcri(CARCRI_SIZE)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics/SIMU_POINT_MAT)
!
! Get list of parameters for integration of constitutive law
!
! --------------------------------------------------------------------------------------------------
!
! Out carcri           : list of parameters for integration of constitutive law
!
! --------------------------------------------------------------------------------------------------
!
    type(Behaviour_PrepCrit) :: behaviourPrepCrit
!
! --------------------------------------------------------------------------------------------------
!
    carcri(1:CARCRI_SIZE) = 0.d0

! - Create carcri informations objects
    call carc_info(behaviourPrepCrit)

! - Read informations from command file
    call carc_read(behaviourPrepCrit)

! - Some checks
    call carc_chck(behaviourPrepCrit)

! - Set in list
    call setBehaviourParaValue(behaviourPrepCrit%v_crit,&
                               behaviourPrepCrit%parm_theta_thm, behaviourPrepCrit%parm_alpha_thm,&
                               behaviourPrepCrit%hho_coef_stab , behaviourPrepCrit%hho_type_stab ,&
                               behaviourPrepCrit%hho_type_calc,&
                               carcriList_ = carcri(1:CARCRI_SIZE))

! - Clean
    deallocate(behaviourPrepCrit%v_crit)
!
end subroutine
