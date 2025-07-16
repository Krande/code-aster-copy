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
#include "asterf_types.h"
!
interface
    subroutine nxacmv(model, materField, mateco, caraElem, listLoad, nume_dof, &
                      solver, l_stat, timeMap, timeParaIn, temp_iter, &
                      vhydr, varc_prev, varc_curr, cn2mbr_stat, &
                      cn2mbr_tran, matass, maprec, cndiri, cncine, &
                      mediri, comporTher, ds_algorom_)
        use Rom_Datastructure_type
        character(len=8), intent(in) :: model, materField, caraElem
        character(len=24), intent(in) :: mateco
        character(len=24), intent(in) :: listLoad
        character(len=24), intent(in) :: nume_dof
        character(len=19), intent(in) :: solver
        character(len=24), intent(in) :: timeMap
        character(len=19), intent(in) :: varc_prev
        character(len=19), intent(in) :: varc_curr
        aster_logical, intent(in) :: l_stat
        real(kind=8), intent(in) :: timeParaIn(6)
        character(len=24), intent(in) :: temp_iter
        character(len=24), intent(in) :: vhydr
        character(len=24), intent(in) :: cn2mbr_stat
        character(len=24), intent(in) :: cn2mbr_tran
        character(len=24), intent(in) :: matass
        character(len=19), intent(in) :: maprec
        character(len=24), intent(in) :: cndiri
        character(len=24), intent(out) :: cncine
        character(len=24), intent(in) :: mediri
        character(len=24), intent(in) :: comporTher
        type(ROM_DS_AlgoPara), optional, intent(in) :: ds_algorom_
    end subroutine nxacmv
end interface
