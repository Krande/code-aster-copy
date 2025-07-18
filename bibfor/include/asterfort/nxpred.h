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
#include "asterf_types.h"
!
interface
    subroutine nxpred(model, mateco, caraElem, listLoad, nume_dof, &
                      solver, l_stat, tpsthe, timeMap, matass, &
                      lonch, maprec, varc_curr, temp_prev, temp_iter, &
                      cn2mbr, hydr_prev, hydr_curr, dry_curr, &
                      comporTher, cndirp, cnchci, vec2nd, vec2ni, &
                      ds_algorom)
        use Rom_Datastructure_type
        character(len=8), intent(in) :: model, caraElem
        character(len=24), intent(in) :: mateco, listLoad
        character(len=24), intent(in) :: nume_dof
        character(len=19), intent(in) :: solver
        real(kind=8), intent(in) :: tpsthe(6)
        character(len=24), intent(in) :: timeMap
        character(len=19), intent(in) :: varc_curr
        aster_logical :: l_stat
        integer(kind=8) :: lonch
        character(len=24) :: matass
        character(len=19) :: maprec
        character(len=24) :: temp_prev
        character(len=24) :: temp_iter
        character(len=24) :: hydr_prev
        character(len=24) :: hydr_curr
        character(len=24) :: dry_curr
        character(len=24) :: comporTher
        character(len=24) :: cndirp
        character(len=24) :: cnchci
        character(len=24), intent(in) :: cn2mbr
        character(len=24) :: vec2nd
        character(len=24) :: vec2ni
        type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    end subroutine nxpred
end interface
