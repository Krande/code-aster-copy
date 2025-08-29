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
subroutine ndxforc_pred(list_func_acti, &
                        model, cara_elem, list_load, &
                        nume_dof, &
                        sddyna, nlDynaDamping, &
                        ds_material, ds_constitutive, ds_system, &
                        ds_measure, sdnume, &
                        sddisc, nume_inst, &
                        hval_incr, hval_algo, &
                        hval_veelem, hval_veasse, &
                        hval_measse, ldccvg)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use NonLinearDyna_module, only: compViteForce
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/ndfdyn.h"
#include "asterfort/diinst.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nonlinDynaMDampCompute.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/nonlinIntForce.h"
#include "asterfort/nonlinLoadCompute.h"
#include "asterfort/nonlinRForceCompute.h"
#include "asterfort/nonlinSubStruCompute.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=24), intent(in) :: model, cara_elem, nume_dof
    character(len=19), intent(in) :: list_load, sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: sddisc, sdnume
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=19), intent(in) :: hval_veelem(*), hval_veasse(*)
    character(len=19), intent(in) :: hval_measse(*)
    integer(kind=8), intent(out) :: ldccvg
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Compute forces for second member at prediction for explicit algorithm
!
! --------------------------------------------------------------------------------------------------
!
! In  list_func_acti   : list of active functionnalities
! In  model            : name of model
! In  cara_elem        : name of elementary characteristics (field)
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  list_load        : name of datastructure for list of loads
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_system        : datastructure for non-linear system management
! IO  ds_measure       : datastructure for measure and statistics management
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current time step
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_veelem      : hat-variable for elementary vectors
! In  hval_veasse      : hat-variable for vectors (node fields)
! In  hval_measse      : hat-variable for matrix
! Out ldccvg           : indicator from integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: phaseType = PRED_EULER
    integer(kind=8), parameter :: iterNewtPred = 0
    integer(kind=8) :: ifm, niv
    character(len=19) :: cndyna, cnsstr, cnhyst
    character(len=19) :: dispCurr, accePrev
    real(kind=8) :: timePrev, timeCurr
    aster_logical :: l_mstp, lDampModal, lDampMatrix, lSuperElement
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_21')
    end if

! - Initializations
    ldccvg = 0

! - Get time
    ASSERT(nume_inst .gt. 0)
    timePrev = diinst(sddisc, nume_inst-1)
    timeCurr = diinst(sddisc, nume_inst)

! - Active functionnalities
    l_mstp = ndynlo(sddyna, 'MULTI_PAS')
    lDampModal = nlDynaDamping%lDampModal
    lDampMatrix = nlDynaDamping%hasMatrDamp
    lSuperElement = isfonc(list_func_acti, 'MACR_ELEM_STAT')

! - Get hat variables
    call nmchex(hval_incr, 'VALINC', 'ACCMOI', accePrev)
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', dispCurr)

! - Compute loads (undead)
    call nonlinLoadCompute('VARI', list_load, &
                           model, cara_elem, nume_dof, list_func_acti, &
                           ds_material, ds_constitutive, ds_measure, &
                           timePrev, timeCurr, &
                           hval_incr, hval_algo, &
                           hval_veelem, hval_veasse)

! - Compute sub-structuring effect on second member
    if (lSuperElement) then
        call nmchex(hval_veasse, 'VEASSE', 'CNSSTR', cnsstr)
        call nonlinSubStruCompute(ds_measure, dispCurr, &
                                  hval_measse, cnsstr)
    end if

! - Compute force for Dirichlet boundary conditions (dualized) - BT.LAMBDA
    call nonlinRForceCompute(model, list_load, &
                             nume_dof, ds_measure, accePrev, &
                             hval_veelem, hval_veasse)

! - Compute effect of dynamic forces (from time discretization scheme)
    call nmchex(hval_veasse, 'VEASSE', 'CNDYNA', cndyna)
    call ndfdyn(sddyna, nlDynaDamping, &
                hval_incr, hval_measse, ds_measure, &
                cndyna)

! - Compute effect of damping (C \cdot \dot{u})
    if (lDampMatrix .and. l_mstp) then
        call nmchex(hval_veasse, 'VEASSE', 'CNHYST', cnhyst)
        call compViteForce(nlDynaDamping, hval_incr, 'VITPLU', cnhyst)
    end if

! - Compute modal damping
    if (lDampModal) then
        call nonlinDynaMDampCompute(phaseType, &
                                    nlDynaDamping, &
                                    nume_dof, ds_measure, &
                                    hval_incr, hval_veasse)
    end if

! - Compute internal forces
    call nonlinIntForce(phaseType, &
                        model, cara_elem, &
                        list_func_acti, iterNewtPred, sdnume, &
                        ds_material, ds_constitutive, &
                        ds_system, ds_measure, &
                        hval_incr, hval_algo, &
                        ldccvg, &
                        sddyna_=sddyna)
!
end subroutine
