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
interface
    subroutine nxrech(model, mateco, caraElem, listLoad, nume_dof, &
                      tpsthe, timeMap, lonch, comporTher, varc_curr, &
                      temp_iter, vtempp, vtempr, temp_prev, hydr_prev, &
                      hydr_curr, vec2nd, cnvabt, &
                      cnresi, rho, iterho, ds_algopara, l_stat)
        use NonLin_Datastructure_type
        character(len=8), intent(in) :: model, caraElem
        character(len=24), intent(in) :: mateco
        character(len=24), intent(in) :: listLoad
        character(len=24), intent(in) :: nume_dof
        real(kind=8), intent(in) :: tpsthe(6)
        character(len=24), intent(in) :: timeMap
        character(len=19), intent(in) :: varc_curr
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        integer(kind=8) :: lonch
        character(len=24) :: comporTher
        character(len=24) :: vtempp
        character(len=24) :: vtempr
        character(len=24) :: temp_prev
        character(len=24) :: temp_iter
        character(len=24) :: hydr_prev
        character(len=24) :: hydr_curr
        character(len=24) :: vec2nd
        character(len=24) :: cnvabt
        character(len=24) :: cnresi
        real(kind=8) :: rho
        integer(kind=8) :: iterho
        aster_logical, intent(in) :: l_stat
    end subroutine nxrech
end interface
