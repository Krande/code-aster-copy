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
subroutine ndasva(sddyna, nlDynaDamping, hval_veasse, cnvady)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nonlinDSVectCombCompute.h"
#include "asterfort/nonlinDSVectCombAddHat.h"
#include "asterfort/nonlinDSVectCombAddDyna.h"
#include "asterfort/nonlinDSVectCombInit.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
!
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19), intent(in) :: hval_veasse(*), cnvady
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Get undead Neumann loads for dynamic
!
! --------------------------------------------------------------------------------------------------
!
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! In  hval_veasse      : hat-variable for vectors (node fields)
! In  cnvady           : name of resultant nodal field
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_mstp, lDampModal, lDampMatrix
    real(kind=8) :: coeam0
    type(NL_DS_VectComb) :: ds_vectcomb
!
! --------------------------------------------------------------------------------------------------
!
    l_mstp = ndynlo(sddyna, 'MULTI_PAS')
    lDampModal = nlDynaDamping%lDampModal
    lDampMatrix = nlDynaDamping%hasMatrDamp

! - Initializations
    call nonlinDSVectCombInit(ds_vectcomb)

! - Coefficients
    coeam0 = ndynre(sddyna, 'COEF_MPAS_FAMO_PREC')

! - Undead dynamic forces
    call nonlinDSVectCombAddHat(hval_veasse, 'CNDYNA', -1.d0, ds_vectcomb)
    if (lDampModal) then
        call nonlinDSVectCombAddHat(hval_veasse, 'CNAMOD', -1.d0, ds_vectcomb)
    end if
    if (lDampMatrix) then
        if (l_mstp) then
            call nonlinDSVectCombAddDyna(sddyna, 'CNHYST', -1.d0*coeam0, ds_vectcomb)
        end if
    end if

! - Combination
    call nonlinDSVectCombCompute(ds_vectcomb, cnvady)
!
end subroutine
