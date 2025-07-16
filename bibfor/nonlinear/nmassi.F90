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
subroutine nmassi(list_func_acti, sddyna, nlDynaDamping, ds_system, hval_incr, hval_veasse, cndonn)
    !
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use NonLinearDyna_module, only: compViteForce
    !
    implicit none
    !
#include "asterf_types.h"
#include "asterfort/infdbg.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmacfi.h"
#include "asterfort/nmacva.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nmchex.h"
#include "asterfort/ndynkk.h"
#include "asterfort/copisd.h"
#include "asterfort/utmess.h"
#include "asterfort/nonlinDSVectCombCompute.h"
#include "asterfort/nonlinDSVectCombAddHat.h"
#include "asterfort/nonlinDSVectCombInit.h"
#include "asterfort/nonlinDSVectCombAddAny.h"
    !
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=19), intent(in) :: sddyna, hval_incr(*), hval_veasse(*)
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: cndonn
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Evaluate second member for initial acceleration
!
! --------------------------------------------------------------------------------------------------
!
! In  list_func_acti   : list of active functionnalities
! In  sddyna           : datastructure for dynamic
! In  nlDynaDamping    : damping parameters
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_veasse      : hat-variable for vectors (nodal fields)
! In  ds_system        : datastructure for non-linear system management
! In  cndonn           : name of nodal field for "given" forces
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: cnffdo, cndfdo, cnfvdo, olhyst, cnhyst
    aster_logical :: l_wave
    aster_logical :: lDampMatrix, lElemDampFromUser
    aster_logical :: l_mstp
    type(NL_DS_VectComb) :: ds_vectcomb
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_17')
    end if
    !
    ! - Active functionnalities
    !
    lDampMatrix = nlDynaDamping%hasMatrDamp
    lElemDampFromUser = nlDynaDamping%lElemDampFromUser
    l_wave = ndynlo(sddyna, 'ONDE_PLANE')
    if (l_wave) then
        call utmess('A', 'MECANONLINE_23')
    end if
    l_mstp = ndynlo(sddyna, 'MULTI_PAS')
    !
    ! - Initializations
    !
    call nonlinDSVectCombInit(ds_vectcomb)
    cnffdo = '&&CNCHAR.FFDO'
    cndfdo = '&&CNCHAR.DFDO'
    cnfvdo = '&&CNCHAR.FVDO'
    !
    ! - Get Dirichlet boundary conditions - B.U
    !
    call nonlinDSVectCombAddHat(hval_veasse, 'CNBUDI', -1.d0, ds_vectcomb)
    !
    ! - Get Neumann and Dirichlet loads
    !
    call nmacfi(list_func_acti, hval_veasse, cnffdo, cndfdo)
    call nonlinDSVectCombAddAny(cndfdo, +1.d0, ds_vectcomb)
    call nonlinDSVectCombAddAny(cnffdo, +1.d0, ds_vectcomb)
    !
    ! - Get variables loads
    !
    call nmacva(hval_veasse, cnfvdo)
    call nonlinDSVectCombAddAny(cnfvdo, +1.d0, ds_vectcomb)
    !
    ! - Compute force induced by damping (C \cdot \dot{u}_{0})
    !
    if (lDampMatrix) then
        if (lElemDampFromUser) then
            call utmess('I', 'MECANONLINE_80')
        else
            call nmchex(hval_veasse, 'VEASSE', 'CNHYST', cnhyst)
            call compViteForce(nlDynaDamping, hval_incr, 'VITMOI', cnhyst)
            call nonlinDSVectCombAddAny(cnhyst, -1.d0, ds_vectcomb)

            ! Save second member for multi-step methods
            if (l_mstp) then
                call ndynkk(sddyna, 'OLDP_CNHYST', olhyst)
                call copisd('CHAMP_GD', 'V', cnhyst, olhyst)
            end if
        end if
    end if
    !
    ! - Add internal forces to second member
    !
    call nonlinDSVectCombAddAny(ds_system%cnfnod, -1.d0, ds_vectcomb)
    !
    ! - Second member (standard)
    !
    call nonlinDSVectCombCompute(ds_vectcomb, cndonn)
    if (niv .ge. 2) then
        call nmdebg('VECT', cndonn, 6)
    end if
    !
end subroutine
