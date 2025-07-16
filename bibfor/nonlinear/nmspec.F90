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
subroutine nmspec(model, ds_material, caraElem, listLoad, listFuncActi, &
                  nume_dof, ds_system, &
                  ds_constitutive, &
                  sddisc, numeTime, &
                  sddyna, sderro, ds_algopara, &
                  ds_measure, &
                  hval_incr, hval_algo, &
                  hval_meelem, &
                  ds_posttimestep)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/isfonc.h"
#include "asterfort/selectListGet.h"
#include "asterfort/nmflam.h"
#include "asterfort/utmess.h"
!
    character(len=24), intent(in) :: model, caraElem
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=19), intent(in) :: listLoad
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: numeTime
    character(len=19), intent(in) :: sddyna
    character(len=24), intent(in) :: sderro
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=19), intent(in) :: hval_meelem(*)
    type(NL_DS_PostTimeStep), intent(inout) :: ds_posttimestep
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Initializations
!
! Spectral analysis
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  ds_material      : datastructure for material parameters
! In  caraElem         : name of elementary characteristics (field)
! In  listLoad        : datastructure for list of loads
! In  listFuncActi     : list of active functionnalities
! In  nume_dof         : name of numbering (NUME_DDL)
! In  ds_constitutive  : datastructure for constitutive laws management
! In  sddisc           : datastructure for time discretization
! In  numeTime        : index of current time step
! In  sddyna           : datastructure for dynamic
! In  sderro           : datastructure for error management (events)
! In  ds_algopara      : datastructure for algorithm parameters
! In  ds_system        : datastructure for non-linear system management
! IO  ds_measure       : datastructure for measure and statistics management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_meelem      : hat-variable for elementary matrix
! IO  ds_posttimestep  : datastructure for post-treatment at each time step
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_mode_vibr, l_crit_stab, l_select
    real(kind=8) :: timeCurr
    character(len=16) :: optionSpec
!
! --------------------------------------------------------------------------------------------------
!
    timeCurr = diinst(sddisc, numeTime)
    l_select = ASTER_FALSE
    optionSpec = ' '

! - Active functionnalites
    l_mode_vibr = isfonc(listFuncActi, 'MODE_VIBR')
    l_crit_stab = isfonc(listFuncActi, 'CRIT_STAB')

! - Compute stability criterion
    if (l_crit_stab) then
        call selectListGet(ds_posttimestep%crit_stab%selector, numeTime, timeCurr, l_select)
        if (l_select) then
            optionSpec = ds_posttimestep%crit_stab%option
! --------- Print
            if (optionSpec .eq. 'FLAMBSTA') then
                call utmess('I', 'MECANONLINE6_2')
            else if (optionSpec .eq. 'FLAMBDYN') then
                call utmess('I', 'MECANONLINE6_2')
            else
                ASSERT(ASTER_FALSE)
            end if
! --------- Compute
            call nmflam(optionSpec, &
                        model, ds_material, caraElem, listLoad, listFuncActi, &
                        nume_dof, ds_system, &
                        ds_constitutive, &
                        sddisc, numeTime, &
                        sddyna, sderro, ds_algopara, &
                        ds_measure, &
                        hval_incr, hval_algo, &
                        hval_meelem, &
                        ds_posttimestep)
        end if
    end if
!
! - Compute vibration modes
!
    if (l_mode_vibr) then
        call selectListGet(ds_posttimestep%mode_vibr%selector, numeTime, timeCurr, l_select)
        if (l_select) then
            optionSpec = ds_posttimestep%mode_vibr%option
! --------- Print
            call utmess('I', 'MECANONLINE6_3')
! --------- Compute
            call nmflam(optionSpec, &
                        model, ds_material, caraElem, listLoad, listFuncActi, &
                        nume_dof, ds_system, &
                        ds_constitutive, &
                        sddisc, numeTime, &
                        sddyna, sderro, ds_algopara, &
                        ds_measure, &
                        hval_incr, hval_algo, &
                        hval_meelem, &
                        ds_posttimestep)
        end if
    end if
!
end subroutine
