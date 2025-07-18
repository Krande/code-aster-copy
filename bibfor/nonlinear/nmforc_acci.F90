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
subroutine nmforc_acci(list_func_acti, &
                       model, cara_elem, nume_dof, &
                       list_load, sddyna, &
                       ds_material, ds_constitutive, ds_system, &
                       ds_measure, ds_inout, &
                       sddisc, nume_inst, &
                       hval_incr, hval_algo, &
                       hval_veelem, hval_veasse, &
                       hval_measse)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/assvec.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nonlinDynaImpeCompute.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/nonlinLoadCompute.h"
#include "asterfort/nonlinLoadDynaCompute.h"
#include "asterfort/nonlinNForceCompute.h"
#include "asterfort/nonlinSubStruCompute.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=24), intent(in) :: model, cara_elem, nume_dof
    character(len=19), intent(in) :: list_load, sddyna
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=19), intent(in) :: hval_veelem(*), hval_veasse(*)
    character(len=19), intent(in) :: hval_measse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Compute forces for second member for initial acceleration
!
! --------------------------------------------------------------------------------------------------
!
! In  list_func_acti   : list of active functionnalities
! In  model            : name of model
! In  cara_elem        : name of elementary characteristics (field)
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  list_load        : name of datastructure for list of loads
! In  sddyna           : name of datastructure for dynamic parameters
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_system        : datastructure for non-linear system management
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_inout         : datastructure for input/output management
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current time step
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_veelem      : hat-variable for elementary vectors
! In  hval_measse      : hat-variable for matrix
! In  hval_veasse      : hat-variable for vectors (node fields)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: phaseTypePred = PRED_EULER
    integer(kind=8) :: ifm, niv
    real(kind=8) :: time_prev, time_curr
    aster_logical :: l_macr, l_impe
    character(len=19) :: disp_curr, cnsstr
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_3')
    end if
!
! - Get time
!
    ASSERT(nume_inst .eq. 0)
    time_prev = 0.d0
    time_curr = diinst(sddisc, nume_inst)
!
! - Active functionnalities
!
    l_impe = ndynlo(sddyna, 'IMPE_ABSO')
    l_macr = isfonc(list_func_acti, 'MACR_ELEM_STAT')
!
! - Compute loads
!
    call nonlinLoadCompute('ACCI', list_load, &
                           model, cara_elem, nume_dof, list_func_acti, &
                           ds_material, ds_constitutive, ds_measure, &
                           time_prev, time_curr, &
                           hval_incr, hval_algo, &
                           hval_veelem, hval_veasse)
!
! - Compute loads (for dynamic)
!
    call nonlinLoadDynaCompute('ACCI', sddyna, &
                               model, nume_dof, &
                               ds_material, ds_measure, ds_inout, &
                               time_prev, time_curr, &
                               hval_veelem, hval_veasse)
!
! - Compute impedance
!
    if (l_impe) then
        call nonlinDynaImpeCompute(phaseTypePred, sddyna, &
                                   model, nume_dof, &
                                   ds_material, ds_measure, &
                                   hval_incr, &
                                   hval_veelem, hval_veasse)
    end if
!
! - Compute sub-structuring effect on second member
!
    if (l_macr) then
        call nmchex(hval_incr, 'VALINC', 'DEPPLU', disp_curr)
        call nmchex(hval_veasse, 'VEASSE', 'CNSSTR', cnsstr)
        call nonlinSubStruCompute(ds_measure, disp_curr, &
                                  hval_measse, cnsstr)
    end if
!
! - Compute nodal force BT . SIGMA (No integration of behaviour)
!
    call nonlinNForceCompute(model, cara_elem, list_func_acti, &
                             ds_material, ds_constitutive, &
                             ds_measure, ds_system, &
                             hval_incr, hval_algo)
    call assvec('V', ds_system%cnfnod, 1, ds_system%vefnod, [1.d0], &
                ds_system%nume_dof)
    if (niv .ge. 2) then
        call nmdebg('VECT', ds_system%cnfnod, 6)
    end if
!
end subroutine
