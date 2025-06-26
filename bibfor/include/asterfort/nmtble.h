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
    subroutine nmtble(loop_exte, model, mesh, ds_contact, &
                      list_func_acti, ds_print, &
                      sderro, ds_conv, sddisc, nume_inst, hval_incr, &
                      hval_algo, ds_algorom)
        use NonLin_Datastructure_type
        use Rom_Datastructure_type
        integer(kind=8), intent(inout) :: loop_exte
        character(len=24), intent(in) :: model
        character(len=8), intent(in) :: mesh
        type(NL_DS_Contact), intent(inout) :: ds_contact
        integer(kind=8), intent(in) :: list_func_acti(*)
        type(NL_DS_Print), intent(inout) :: ds_print
        character(len=24), intent(in) :: sderro
        type(NL_DS_Conv), intent(in) :: ds_conv
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: nume_inst
        character(len=19), intent(in) :: hval_incr(*)
        character(len=19), intent(in) :: hval_algo(*)
        type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
    end subroutine nmtble
end interface
