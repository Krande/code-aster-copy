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

subroutine cfmmcv(mesh, model_, list_func_acti, iter_newt, nume_inst, &
                  sddyna, ds_measure, sddisc, sderro, hval_incr, &
                  hval_algo, ds_print, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfconv.h"
#include "asterfort/isfonc.h"
#include "asterfort/mm_cycl_print.h"
#include "asterfort/mmbclc.h"
#include "asterfort/mmbouc.h"
#include "asterfort/nmcrel.h"
#include "asterfort/nmimci.h"
#include "asterfort/nmimck.h"
#include "asterfort/nmimcr.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: model_
    integer(kind=8), intent(in) :: list_func_acti(*)
    integer(kind=8), intent(in) :: iter_newt
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: sddisc
    character(len=19), intent(in) :: sddyna
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=24), intent(in) :: sderro
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: hval_algo(*)
    type(NL_DS_Print), intent(inout) :: ds_print
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! All methods - Evaluate convergence
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  iter_newt        : index of current Newton iteration
! In  nume_inst        : index of current time step
! In  sddisc           : datastructure for time discretization
! In  sddyna           : dynamic parameters datastructure
! IO  ds_measure       : datastructure for measure and statistics management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  sderro           : datastructure for errors during algorithm
! In  hval_algo        : hat-variable for algorithms fields
! IO  ds_print         : datastructure for printing parameters
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_cont_disc, l_cont_cont, l_newt_cont
    aster_logical :: loop_cont_conv, loop_geom_error, l_all_verif
    character(len=8) :: model
    real(kind=8) :: r8bid, loop_cont_vale
    integer(kind=8) :: loop_cont_vali
!
! --------------------------------------------------------------------------------------------------
!
    r8bid = 0.d0
    model = model_(1:8)

! - Get parameters
    l_cont_disc = isfonc(list_func_acti, 'CONT_DISCRET')
    l_cont_cont = isfonc(list_func_acti, 'CONT_CONTINU')
    l_newt_cont = isfonc(list_func_acti, 'CONT_NEWTON')
    l_all_verif = isfonc(list_func_acti, 'CONT_ALL_VERIF')

! - Values in convergence table: not affected
    call nmimck(ds_print, 'BOUC_NOEU', ' ', .false._1)
    call nmimcr(ds_print, 'BOUC_VALE', r8bid, .false._1)

! - Convergence for contact discrete methods
    if (l_cont_disc) then
        call cfconv(mesh, ds_measure, sderro, hval_algo, ds_print, &
                    ds_contact)
    end if
!
! - Applying generalized Newton method at Newton's iteration
!
    if (l_newt_cont) then
        call mmbclc(mesh, model, iter_newt, nume_inst, ds_measure, &
                    sddisc, sddyna, hval_incr, hval_algo, ds_contact)

        call mmbouc(ds_contact, 'Cont', 'Is_Convergence', loop_state_=loop_cont_conv)

        if (loop_cont_conv) then
            call nmcrel(sderro, 'DIVE_CTCC', .false._1)
        else
            call nmcrel(sderro, 'DIVE_CTCC', .true._1)
        end if

        call mmbouc(ds_contact, 'Geom', 'Is_Error', loop_state_=loop_geom_error)
        if (loop_geom_error) then
            call nmcrel(sderro, 'ERRE_APPA', .true._1)
        end if

        if (.not. l_all_verif) then
            call mmbouc(ds_contact, 'Cont', 'Get_Vale', loop_vale_=loop_cont_vale)
            loop_cont_vali = nint(loop_cont_vale)
            call nmimci(ds_print, 'CONT_NEWT', loop_cont_vali, .true._1)
        end if
    end if
!
! - Cycling informations printing in convergence table
!
    if (l_cont_cont) then
        ds_print%resi_pressure = ds_contact%resi_press_glob
        call mm_cycl_print(ds_print, ds_measure)
    end if
!
end subroutine
