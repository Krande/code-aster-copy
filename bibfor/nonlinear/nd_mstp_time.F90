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

subroutine nd_mstp_time(ds_inout, list_func_acti, time_prev_step, l_comp_mstp)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/isfonc.h"
#include "asterc/r8vide.h"
#include "asterfort/utmess.h"
#include "asterfort/rs_getlast.h"
#include "asterfort/rsadpa.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_InOut), intent(in) :: ds_inout
    integer(kind=8), intent(in) :: list_func_acti(*)
    real(kind=8), intent(out) :: time_prev_step
    aster_logical, intent(out) :: l_comp_mstp
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Dynamic
!
! Get previous time for multi-step schemes
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_inout         : datastructure for input/output management
! In  list_func_acti   : list of active functionnalities
! Out time_prev_step   : previous time for multi-step schemes
! Out l_comp_mstp      : .true. if compute second member
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_reuse, l_init_state, l_stin_evol
    integer(kind=8) :: init_nume, nume_prev_step, nume_last
    integer(kind=8) :: jv_para
    character(len=24) :: stin_evol
    character(len=8) :: result
!
! --------------------------------------------------------------------------------------------------
!
    time_prev_step = r8vide()
    l_comp_mstp = .false.
    result = ds_inout%result
    init_nume = ds_inout%init_nume
!
! - Does ETAT_INIT (initial state) exist ?
!
    l_init_state = isfonc(list_func_acti, 'ETAT_INIT')
!
! - Reuse previous results ?
!
    l_reuse = isfonc(list_func_acti, 'REUSE')
!
! - Get name of result datastructure in ETAT_INIT
!
    l_stin_evol = ds_inout%l_stin_evol
    stin_evol = ds_inout%stin_evol
!
! - Initial state: get time if possible
!
    if (l_init_state) then
!
        nume_prev_step = init_nume-1
!
! ----- Get previous time
!
        if (l_stin_evol) then
            if (nume_prev_step .le. 0) then
                call utmess('I', 'DYNAMIQUE_50')
            else
                call rsadpa(stin_evol, 'L', 1, 'INST_PREC', nume_prev_step, &
                            0, sjv=jv_para, istop=0)
                time_prev_step = zr(jv_para)
                if (time_prev_step .eq. r8vide()) then
                    call utmess('I', 'DYNAMIQUE_51')
                else
                    l_comp_mstp = .true.
                end if
            end if
        else
            call utmess('I', 'DYNAMIQUE_53')
        end if
    end if
!
! - Reuse old results: get time if possible
!
    if (l_reuse) then
        call rs_getlast(result, nume_last)
        call rsadpa(result, 'L', 1, 'INST_PREC', nume_last, &
                    0, sjv=jv_para, istop=0)
        time_prev_step = zr(jv_para)
        if (time_prev_step .eq. r8vide()) then
            call utmess('I', 'DYNAMIQUE_51')
        else
            l_comp_mstp = .true.
        end if
    end if
!
end subroutine
