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
! aslint: disable=W1403
!
subroutine nonlinDSMaterialCreate(ds_material)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
!
    type(NL_DS_Material), intent(out) :: ds_material
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Constitutive laws
!
! Create material management datastructure
!
! --------------------------------------------------------------------------------------------------
!
! Out ds_material      : datastructure for material parameters
!
! --------------------------------------------------------------------------------------------------
!
    ds_material%mater = ' '
    ds_material%mateco = ' '
    ds_material%varc_refe = '&&OP0070.VARC_REFE'
    ds_material%fvarc_init = '&&OP0070.FVARC_INIT'
    ds_material%fvarc_pred = '&&OP0070.FVARC_PRED'
    ds_material%fvarc_curr = '&&OP0070.FVARC_CURR'
!
end subroutine
