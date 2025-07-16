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
    subroutine accel0(model, numeDof, listFuncActi, listLoad,&
                      ds_contact, maprec, solveu, &
                      sddyna, nlDynaDamping,&
                      ds_measure, ds_system,&
                      hval_meelem, hval_measse,&
                      hval_veelem, hval_veasse,&
                      hval_incr, hval_algo)
        use NonLin_Datastructure_type
        use NonLinearDyna_type
        character(len=19), intent(in) :: solveu, maprec, listLoad
        character(len=19), intent(in) :: sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        type(NL_DS_Measure), intent(inout) :: ds_measure
        character(len=24), intent(in) :: numeDof, model
        type(NL_DS_Contact), intent(in) :: ds_contact
        character(len=19), intent(in) :: hval_meelem(*), hval_measse(*)
        character(len=19), intent(in) :: hval_veasse(*), hval_veelem(*)
        character(len=19), intent(in) :: hval_algo(*), hval_incr(*)
        integer(kind=8), intent(in) :: listFuncActi(*)
        type(NL_DS_System), intent(in) :: ds_system
    end subroutine accel0
end interface
