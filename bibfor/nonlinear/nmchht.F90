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
subroutine nmchht(model, ds_material, cara_elem, ds_constitutive, &
                  list_load, nume_dof, list_func_acti, ds_measure, &
                  sddyna, sddisc, sdnume, &
                  hval_incr, hval_algo, hval_measse, ds_inout)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/diinst.h"
#include "asterfort/isfonc.h"
#include "asterfort/mecact.h"
#include "asterfort/nd_mstp_time.h"
#include "asterfort/ndynkk.h"
#include "asterfort/nmcha0.h"
#include "asterfort/nmchai.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdidi.h"
#include "asterfort/nmvcex.h"
#include "asterfort/nmvcle.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/nonlinIntForce.h"
#include "asterfort/nonlinLoadCompute.h"
#include "asterfort/nonlinLoadDynaCompute.h"
#include "asterfort/nonlinSubStruCompute.h"
#include "asterfort/utmess.h"
!
    character(len=24), intent(in) :: model
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=24), intent(in) :: cara_elem
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=24), intent(in) :: nume_dof
    character(len=19), intent(in) :: list_load
    integer(kind=8), intent(in) :: list_func_acti(*)
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: sddyna, sddisc, sdnume
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*), hval_measse(*)
    type(NL_DS_InOut), intent(in) :: ds_inout
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Dynamic
!
! Compute previous second member for multi-step schemes
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  ds_material      : datastructure for material parameters
! In  cara_elem        : name of elementary characteristics (field)
! In  ds_constitutive  : datastructure for constitutive laws management
! In  nume_dof         : name of numbering (NUME_DDL)
! In  list_load        : name of datastructure for list of loads
! In  list_func_acti   : list of active functionnalities
! IO  ds_measure       : datastructure for measure and statistics management
! In  sddyna           : dynamic parameters datastructure
! In  sddisc           : datastructure for time discretization
! In  sdnume           : datastructure for dof positions
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_measse      : hat-variable for matrix
! In  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: phaseType = CORR_NEWTON
    integer(kind=8), parameter:: zveass = 19, zveelm = 12
    character(len=19) :: hval_veelem(zveelm)
    character(len=19) :: hval_veasse(zveass)
    aster_logical :: l_didi, l_comp_mstp, l_macr
    character(len=19) :: vefint, vedido
    character(len=19) :: vefedo, veondp, vedidi, vesstf
    character(len=19) :: cnfedo, cndidi, cnfint
    character(len=19) :: cndido, cncine, cnviss, cnhyst
    character(len=19) :: cnondp, cnsstf, cnsstr
    character(len=19) :: disp_prev
    character(len=19) :: varc_prev, varc_curr, time_prev, time_curr
    real(kind=8) :: time_init, time_prev_step
    integer(kind=8) :: iter_newt, ldccvg, nmax
    character(len=4) :: mode
    type(NL_DS_System) :: ds_system
!
! --------------------------------------------------------------------------------------------------
!
    iter_newt = 0
    call nmchai('VEELEM', 'LONMAX', nmax)
    ASSERT(nmax .eq. zveelm)
    call nmchai('VEASSE', 'LONMAX', nmax)
    ASSERT(nmax .eq. zveass)

! - Active functionnalities
    l_didi = isfonc(list_func_acti, 'DIDI')
    l_macr = isfonc(list_func_acti, 'MACR_ELEM_STAT')

! - Initial time
    time_init = diinst(sddisc, 0)

! - Get previous time
    call nd_mstp_time(ds_inout, list_func_acti, time_prev_step, l_comp_mstp)

! - Protection
    if (abs(time_prev_step-time_init) .le. r8prem()) then
        l_comp_mstp = ASTER_FALSE
        call utmess('A', 'DYNAMIQUE_52')
    end if

! - Computation
    if (l_comp_mstp) then
! ----- Check Delta_time
        if (time_init-time_prev_step .lt. 0.d0) then
            call utmess('F', 'DYNAMIQUE1_6')
        end if

! ----- Create <CARTE> for time
        call nmchex(hval_incr, 'VALINC', 'COMMOI', varc_prev)
        call nmchex(hval_incr, 'VALINC', 'COMPLU', varc_curr)
        call nmvcex('INST', varc_prev, time_prev)
        call copisd('CHAMP_GD', 'V', time_prev, varc_curr(1:14)//'.INST')
        call nmvcex('INST', varc_curr, time_curr)
        call mecact('V', time_prev, 'MODELE', model(1:8)//'.MODELE', 'INST_R', &
                    ncmp=1, nomcmp='INST', sr=time_prev_step)
        call mecact('V', time_curr, 'MODELE', model(1:8)//'.MODELE', 'INST_R', &
                    ncmp=1, nomcmp='INST', sr=time_init)
! ----- Create external state variables
        call nmvcle(model, ds_material%mater, cara_elem, time_prev_step, varc_prev)
        call nmvcle(model, ds_material%mater, cara_elem, time_init, varc_curr)

! ----- Get fields from hat-variables -
        call ndynkk(sddyna, 'OLDP_VEFEDO', vefedo)
        call ndynkk(sddyna, 'OLDP_VEDIDO', vedido)
        call ndynkk(sddyna, 'OLDP_VEDIDI', vedidi)
        call ndynkk(sddyna, 'OLDP_VEFINT', vefint)
        call ndynkk(sddyna, 'OLDP_VEONDP', veondp)
        call ndynkk(sddyna, 'OLDP_VESSTF', vesstf)
        call ndynkk(sddyna, 'OLDP_CNFEDO', cnfedo)
        call ndynkk(sddyna, 'OLDP_CNDIDO', cndido)
        call ndynkk(sddyna, 'OLDP_CNDIDI', cndidi)
        call ndynkk(sddyna, 'OLDP_CNFINT', cnfint)
        call ndynkk(sddyna, 'OLDP_CNONDP', cnondp)
        call ndynkk(sddyna, 'OLDP_CNCINE', cncine)
        call ndynkk(sddyna, 'OLDP_CNVISS', cnviss)
        call ndynkk(sddyna, 'OLDP_CNHYST', cnhyst)
        call ndynkk(sddyna, 'OLDP_CNSSTF', cnsstf)
        call ndynkk(sddyna, 'OLDP_CNSSTR', cnsstr)
! ----- New hat variables
        call nmcha0('VEELEM', 'ALLINI', ' ', hval_veelem)
        call nmcha0('VEELEM', 'CNFEDO', vefedo, hval_veelem)
        call nmcha0('VEELEM', 'CNDIDO', vedido, hval_veelem)
        call nmcha0('VEELEM', 'CNDIDI', vedidi, hval_veelem)
        call nmcha0('VEELEM', 'CNONDP', veondp, hval_veelem)
        call nmcha0('VEELEM', 'CNSSTF', vesstf, hval_veelem)
        call nmcha0('VEASSE', 'ALLINI', ' ', hval_veasse)
        call nmcha0('VEASSE', 'CNFEDO', cnfedo, hval_veasse)
        call nmcha0('VEASSE', 'CNDIDO', cndido, hval_veasse)
        call nmcha0('VEASSE', 'CNDIDI', cndidi, hval_veasse)
        call nmcha0('VEASSE', 'CNONDP', cnondp, hval_veasse)
        call nmcha0('VEASSE', 'CNCINE', cncine, hval_veasse)
        call nmcha0('VEASSE', 'CNVISS', cnviss, hval_veasse)
        call nmcha0('VEASSE', 'CNHYST', cnhyst, hval_veasse)
        call nmcha0('VEASSE', 'CNSSTF', cnsstf, hval_veasse)
        call nmcha0('VEASSE', 'CNSSTR', cnsstr, hval_veasse)

! ----- Set in the datastructure ds_system
        ds_system%veinte = vefint
        ds_system%cnfint = cnfint
        ds_system%nume_dof = nume_dof
! ----- Compute forces from macro-elements
        if (l_macr) then
            call nmchex(hval_incr, 'VALINC', 'DEPMOI', disp_prev)
            call nonlinSubStruCompute(ds_measure, disp_prev, &
                                      hval_measse, cnsstr)
        end if
! ----- Compute internal forces
        call nonlinIntForce(phaseType, &
                            model, cara_elem, &
                            list_func_acti, iter_newt, sdnume, &
                            ds_material, ds_constitutive, &
                            ds_system, ds_measure, &
                            hval_incr, hval_algo, &
                            ldccvg)
! ----- Compute forces
        mode = 'FIXE'
        call nonlinLoadCompute(mode, list_load, &
                               model, cara_elem, nume_dof, list_func_acti, &
                               ds_material, ds_constitutive, ds_measure, &
                               time_prev_step, time_init, &
                               hval_incr, hval_algo, &
                               hval_veelem, hval_veasse)
        call nonlinLoadDynaCompute(mode, sddyna, &
                                   model, nume_dof, &
                                   ds_material, ds_measure, ds_inout, &
                                   time_prev_step, time_init, &
                                   hval_veelem, hval_veasse)
! ----- Compute vector for DIDI loads
        if (l_didi) then
            call nmdidi(ds_inout, model, list_load, nume_dof, hval_incr, &
                        hval_veelem, hval_veasse)
        end if
    end if
end subroutine
