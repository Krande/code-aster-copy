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

subroutine nmaffi(list_func_acti, ds_conv, ds_print, sderro, sddisc, &
                  loop_name)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmaffm.h"
#include "asterfort/nmerim.h"
#include "asterfort/nmevim.h"
#include "asterfort/nmimpr.h"
#include "asterfort/nmimps.h"
#include "asterfort/nmimpx.h"
#include "asterfort/nmlecv.h"
#include "asterfort/nmltev.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer(kind=8), intent(in) :: list_func_acti(*)
    type(NL_DS_Conv), intent(in) :: ds_conv
    type(NL_DS_Print), intent(inout) :: ds_print
    character(len=24), intent(in) :: sderro
    character(len=19), intent(in) :: sddisc
    character(len=4), intent(in) :: loop_name
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Print management
!
! Print during loop
!
! --------------------------------------------------------------------------------------------------
!
! In  list_func_acti   : list of active functionnalities
! IO  ds_print         : datastructure for printing parameters
! In  ds_conv          : datastructure for convergence management
! In  sderro           : name of datastructure for error management (events)
! In  sddisc           : name of datastructure for time discretization
! In  loop_name        : name of loop
!                         'NEWT' - Newton loop
!                         'FIXE' - Fixed points loop
!                         'INST' - Step time loop
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_error
    aster_logical :: cvnewt, cvinst
    aster_logical :: l_line_print
    aster_logical :: l_loop_cont, l_dyna_expl
!
! --------------------------------------------------------------------------------------------------
!
    l_line_print = .false.
!
! - Active functionnalites
!
    l_loop_cont = isfonc(list_func_acti, 'BOUCLE_EXTERNE')
    l_dyna_expl = isfonc(list_func_acti, 'EXPLICITE')
!
! - Convergence state of loops
!
    call nmlecv(sderro, 'NEWT', cvnewt)
    call nmlecv(sderro, 'INST', cvinst)
!
! - Set marks in rows
!
    call nmaffm(sderro, ds_print, loop_name)
!
! - Is error event occurred ?
!
    call nmltev(sderro, 'ERRI', loop_name, l_error)
!
! - Print line of convergence table ?
!
    if (loop_name .eq. 'NEWT') then
        if (cvnewt) then
            if (l_loop_cont) then
                l_line_print = .false.
            else
                l_line_print = .true.
            end if
        else
            l_line_print = .true.
        end if
    else if (loop_name .eq. 'FIXE') then
        if (l_loop_cont) then
            l_line_print = .true.
        end if
        if (.not. cvnewt) then
            l_line_print = .false.
        end if
    else if (loop_name .eq. 'INST') then
        l_line_print = .false.
    end if
!
! - Print line in convergence table
!
    if (l_line_print) then
        call nmimpr(ds_print)
    end if
!
! - Print separator line in convergence table
!
    if (l_line_print) then
        if (cvnewt .and. .not. (l_error)) then
            if (ds_print%l_print) then
                call nmimpx(ds_print)
            end if
        end if
    end if
!
! - Print error
!
    if (l_error) then
        call nmimpx(ds_print)
        call nmerim(sderro)
    end if
!
! - Print event messages
!
    call nmevim(ds_print, sddisc, sderro, loop_name)
!
! - Print residuals summary at end of step
!
    if (cvinst .and. .not. l_dyna_expl) then
        call nmimps(ds_print, ds_conv, sderro)
    end if
!
end subroutine
