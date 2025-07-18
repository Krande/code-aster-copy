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
subroutine nmimin(list_func_acti, sddisc, sdsuiv, nume_inst, ds_print)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/InitTableCvg.h"
#include "asterfort/nonlinDSPrintInitTimeStep.h"
#include "asterfort/nonlinDSPrintHeadTimeStep.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: list_func_acti(*)
    character(len=19), intent(in) :: sddisc
    character(len=24), intent(in) :: sdsuiv
    integer(kind=8), intent(in) :: nume_inst
    type(NL_DS_Print), intent(inout) :: ds_print
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Print management
!
! Initializations for new step time
!
! --------------------------------------------------------------------------------------------------
!
! In  list_func_acti   : list of active functionnalities
! In  sddisc           : name of datastructure for time discretization
! In  sdsuiv           : datastructure for DOF monitoring
! In  nume_inst        : index of current time step
! IO  ds_print         : datastructure for printing parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_31')
    end if
!
! - Activations for convergence table
!
    call InitTableCvg(list_func_acti, sdsuiv, ds_print)
!
! - Initializations for convergence table
!
    call nonlinDSPrintInitTimeStep(ds_print)
!
! - Print head for new step time
!
    call nonlinDSPrintHeadTimeStep(sddisc, nume_inst, ds_print)
!
end subroutine
