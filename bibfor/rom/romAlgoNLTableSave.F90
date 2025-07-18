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
subroutine romAlgoNLTableSave(nume_store, time_curr, paraAlgo)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jeveuo.h"
#include "asterfort/romTableSave.h"
!
    integer(kind=8), intent(in) :: nume_store
    real(kind=8), intent(in) :: time_curr
    type(ROM_DS_AlgoPara), intent(in) :: paraAlgo
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Solving non-linear problem
!
! Save table for the reduced coordinates
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_store       : index to store in results
! In  time_curr        : current time
! In  paraAlgo         : datastructure for ROM parameters
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: gamma
    integer(kind=8) :: nbMode
    real(kind=8), pointer :: v_gamma(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    gamma = paraAlgo%gamma
    nbMode = paraAlgo%ds_empi%nbMode
!
! - Access to reduced coordinates
!
    call jeveuo(gamma, 'L', vr=v_gamma)
!
! - Save in table
!
    call romTableSave(paraAlgo%tablResu, nbMode, v_gamma, &
                      nume_store, time_curr)
!
end subroutine
