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
! aslint: disable=W1403
!
subroutine romGreedyAlgoDSInit(ds_solveDOM, ds_solveROM, ds_algoGreedy)
!
    use Rom_Datastructure_type
!
    implicit none
!
    type(ROM_DS_Solve), intent(in)       :: ds_solveDOM
    type(ROM_DS_Solve), intent(in)       :: ds_solveROM
    type(ROM_DS_AlgoGreedy), intent(out) :: ds_algoGreedy
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Initializations
!
! Initialisation of datastructure for greedy algorithm
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_solveDOM      : datastructure for datastructure to solve systems (DOM)
! In  ds_solveROM      : datastructure for datastructure to solve systems (ROM)
! Out ds_algoGreedy    : datastructure for Greedy algorithm
!
! --------------------------------------------------------------------------------------------------
!
    ds_algoGreedy%coef_redu = '&&OP0053.COEF_REDU'
    ds_algoGreedy%resi_vect = '&&OP0053.RESI_VECT'
    ds_algoGreedy%solveDOM = ds_solveDOM
    ds_algoGreedy%solveROM = ds_solveROM
!
end subroutine
