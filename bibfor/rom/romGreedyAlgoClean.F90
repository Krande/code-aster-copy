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
! aslint: disable=W1403
! person_in_charge: mickael.abbas at edf.fr
!
subroutine romGreedyAlgoClean(ds_algoGreedy)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterfort/as_deallocate.h"
!
    type(ROM_DS_AlgoGreedy), intent(in) :: ds_algoGreedy
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction
!
! Clean datastructure for greedy algorithm
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_algoGreedy    : datastructure for Greedy algorithm
!
! --------------------------------------------------------------------------------------------------
!
    AS_DEALLOCATE(vr=ds_algoGreedy%resi_norm)
!
end subroutine
