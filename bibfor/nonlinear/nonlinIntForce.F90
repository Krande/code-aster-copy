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
subroutine nonlinIntForce(phaseType, &
                          model, cara_elem, &
                          list_func_acti, iter_newt, sdnume, &
                          ds_material, ds_constitutive, &
                          ds_system, ds_measure, &
                          hval_incr, hval_algo, &
                          ldccvg, &
                          sddyna_, &
                          ds_algorom_)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
    use NonLinear_module, only: inteForceGetOption
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/assert.h"
#include "asterfort/nmfint.h"
#include "asterfort/nonlinNForceCompute.h"
#include "asterfort/nonlinIntForceAsse.h"
!
    integer(kind=8), intent(in) :: phaseType
    character(len=24), intent(in) :: model, cara_elem
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=19), intent(in) :: sdnume
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Measure), intent(inout) :: ds_measure
    integer(kind=8), intent(in) :: iter_newt
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    integer(kind=8), intent(out) :: ldccvg
    character(len=19), optional, intent(in) :: sddyna_
    type(ROM_DS_AlgoPara), optional, intent(in) :: ds_algorom_
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Compute internal forces
!
! --------------------------------------------------------------------------------------------------
!
! In  phaseType        : name of current phase (prediction/correction/internal forces)
! In  model            : name of model
! In  cara_elem        : name of elementary characteristics (field)
! In  list_func_acti   : list of active functionnalities
! In  sddyna           : datastructure for dynamic
! In  sdnume           : datastructure for dof positions
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_system        : datastructure for non-linear system management
! IO  ds_measure       : datastructure for measure and statistics management
! In  time_prev        : time at beginning of time step
! In  time_curr        : time at end of time step
! In  iter_newt        : index of current Newton iteration
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  ds_algorom       : datastructure for ROM parameters
! Out ldccvg           : indicator from integration of behaviour
!                -1 : PAS D'INTEGRATION DU COMPORTEMENT
!                 0 : CAS DE FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE L'INTEGRATION DE LA LDC
!                 3 : SIZZ PAS NUL POUR C_PLAN DEBORST
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: sddyna
    type(ROM_DS_AlgoPara) :: ds_algorom
    aster_logical :: lNodeComp, lInteComp
    integer(kind=8) :: typeAsse
!
! --------------------------------------------------------------------------------------------------
!
    sddyna = ' '
    if (present(sddyna_)) then
        sddyna = sddyna_
    end if
    if (present(ds_algorom_)) then
        ds_algorom = ds_algorom_
    end if
    if (phaseType .eq. PRED_EULER) then
        ASSERT(iter_newt .eq. 0)
    end if
!
! - Get options to compute internal forces
!
    call inteForceGetOption(phaseType, list_func_acti, ds_algorom, &
                            lNodeComp, lInteComp, typeAsse)
!
! - Direct computation (no integration of behaviour)
!
    if (lNodeComp) then
        call nonlinNForceCompute(model, cara_elem, list_func_acti, &
                                 ds_material, ds_constitutive, &
                                 ds_measure, ds_system, &
                                 hval_incr, hval_algo)
        ldccvg = -1
    end if
!
! - Integration of behaviour
!
    if (lInteComp) then
        call nmfint(model, cara_elem, &
                    ds_material, ds_constitutive, &
                    list_func_acti, iter_newt, ds_measure, ds_system, &
                    hval_incr, hval_algo, &
                    ldccvg, sddyna)
    end if
!
! - Assembly
!
    if (ldccvg .ne. 1) then
        call nonlinIntForceAsse(typeAsse, list_func_acti, sdnume, &
                                ds_material, ds_constitutive, ds_system)
    end if
!
end subroutine
