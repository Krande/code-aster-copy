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
    subroutine nxinit(mesh, model, materField, &
                      caraElem, compor, listLoad, &
                      para, nume_dof, &
                      sddisc, ds_inout, sdobse, &
                      time, ds_algopara, &
                      ds_algorom, ds_print, vhydr, &
                      l_stat, l_evol, l_rom, &
                      l_line_search, lnkry, l_dry)
        use NonLin_Datastructure_type
        use Rom_Datastructure_type
        character(len=24), intent(in) ::  compor
        character(len=24), intent(in) :: listLoad
        real(kind=8), intent(in) :: para(*)
        character(len=24), intent(out) :: nume_dof
        character(len=8), intent(in) :: materField, caraElem, model, mesh
        character(len=19), intent(in) :: sddisc
        type(NL_DS_InOut), intent(inout) :: ds_inout
        character(len=19), intent(out) :: sdobse
        character(len=24), intent(out) :: time
        type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
        type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
        type(NL_DS_Print), intent(inout) :: ds_print
        character(len=24), intent(in) :: vhydr
        aster_logical, intent(out) :: l_stat, l_evol, l_rom, l_line_search, lnkry
        aster_logical, intent(in) :: l_dry
    end subroutine nxinit
end interface
