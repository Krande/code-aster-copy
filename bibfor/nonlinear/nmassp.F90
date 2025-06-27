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
subroutine nmassp(listFuncActi, &
                  sddyna, nlDynaDamping, &
                  ds_system, ds_contact, hval_veasse, &
                  cnpilo, cndonn)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/ndassp.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nsassp.h"
#include "asterfort/vtzero.h"
!
    integer(kind=8), intent(in) :: listFuncActi(*)
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: hval_veasse(*)
    character(len=19), intent(in) :: cnpilo, cndonn
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Evaluate second member for prediction
!
! --------------------------------------------------------------------------------------------------
!
! In  listFuncActi     : list of active functionnalities
! In  ds_system        : datastructure for non-linear system management
! In  ds_contact       : datastructure for contact management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! In  hval_veasse      : hat-variable for vectors (node fields)
! In  cndonn           : name of nodal field for "given" forces
! In  cnpilo           : name of nodal field for "pilotage" forces
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_stat, l_dyna
!
! --------------------------------------------------------------------------------------------------
!
    call vtzero(cnpilo)
    call vtzero(cndonn)
!
! - Active functionnalities
!
    l_stat = ndynlo(sddyna, 'STATIQUE')
    l_dyna = ndynlo(sddyna, 'DYNAMIQUE')
!
! - Evaluate second member for prediction
!
    if (l_dyna) then
        call ndassp(listFuncActi, ds_contact, ds_system, &
                    sddyna, nlDynaDamping, &
                    hval_veasse, cndonn)
    else if (l_stat) then
        call nsassp(listFuncActi, ds_contact, ds_system, &
                    hval_veasse, cnpilo, cndonn)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
