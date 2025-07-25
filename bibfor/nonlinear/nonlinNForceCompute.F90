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
subroutine nonlinNForceCompute(model, cara_elem, list_func_acti, &
                               ds_material, ds_constitutive, &
                               ds_measure, ds_system, &
                               hval_incr, hval_algo)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmtime.h"
#include "asterfort/nmvcex.h"
#include "asterfort/utmess.h"
#include "asterfort/vefnme.h"
!
    character(len=24), intent(in) :: model, cara_elem
    integer(kind=8), intent(in) :: list_func_acti(*)
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Compute nodal force BT . SIGMA (No integration of behaviour)
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  cara_elem        : name of elementary characteristics (field)
! In  list_func_acti   : list of active functionnalities
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_system        : datastructure for non-linear system management
! In  time_prev        : time at beginning of time step
! In  time_curr        : time at end of time step
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: disp_prev, strx_prev, sigm_prev, varc_prev
    character(len=19) :: disp_cumu_inst, sigm_extr
    character(len=24) :: vrcmoi
    character(len=16), parameter :: option = 'FORC_NODA'
    aster_logical :: l_implex
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_8')
    end if

! - Hat variable
    call nmchex(hval_algo, 'SOLALG', 'DEPDEL', disp_cumu_inst)
    call nmchex(hval_incr, 'VALINC', 'DEPMOI', disp_prev)
    call nmchex(hval_incr, 'VALINC', 'STRMOI', strx_prev)
    call nmchex(hval_incr, 'VALINC', 'SIGMOI', sigm_prev)
    call nmchex(hval_incr, 'VALINC', 'SIGEXT', sigm_extr)
    call nmchex(hval_incr, 'VALINC', 'COMMOI', varc_prev)
    call nmvcex('TOUT', varc_prev, vrcmoi)

! - Active functionnalities
    l_implex = isfonc(list_func_acti, 'IMPLEX')

! - Launch timer
    call nmtime(ds_measure, 'Init', '2nd_Member')
    call nmtime(ds_measure, 'Launch', '2nd_Member')

! - Compute
    if (l_implex) then
        call vefnme(option, model, ds_material%mateco, cara_elem, &
                    ds_constitutive%compor, 0, ' ', &
                    vrcmoi, sigm_extr, ' ', &
                    disp_prev, &
                    'V', ds_system%vefnod)
    else
        call vefnme(option, model, ds_material%mateco, cara_elem, &
                    ds_constitutive%compor, 0, ' ', &
                    vrcmoi, sigm_prev, strx_prev, &
                    disp_prev, &
                    'V', ds_system%vefnod)
    end if

! - Stop timer
    call nmtime(ds_measure, 'Stop', '2nd_Member')
!
end subroutine
