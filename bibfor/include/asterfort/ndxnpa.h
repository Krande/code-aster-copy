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
    subroutine ndxnpa(model         , cara_elem,&
                      list_func_acti, ds_print,&
                      ds_material   , ds_constitutive,&
                      sddisc, sddyna, sdnume, nume_dof, nume_inst  ,&
                      hval_incr     , hval_algo)
        use NonLin_Datastructure_type
        character(len=24), intent(in) :: model, cara_elem
        integer(kind=8), intent(in) :: list_func_acti(*)
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        character(len=19), intent(in) :: sddyna, sdnume, sddisc
        type(NL_DS_Print), intent(inout) :: ds_print
        integer(kind=8), intent(in) :: nume_inst
        character(len=24), intent(in) :: nume_dof
        character(len=19), intent(in) :: hval_algo(*), hval_incr(*)
    end subroutine ndxnpa
end interface
