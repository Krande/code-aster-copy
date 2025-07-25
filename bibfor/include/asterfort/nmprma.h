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
    subroutine nmprma(listFuncActi, &
                      modelz, caraElem, &
                      ds_material, ds_constitutive, &
                      listLoad, sddyna, nlDynaDamping, &
                      sddisc, numeTime, &
                      ds_algopara, ds_contact, ds_algorom, &
                      ds_print, ds_measure, &
                      hval_incr, hval_algo, &
                      hval_meelem, hval_measse, &
                      numeDof, numeDofFixe, &
                      solveu, ds_system, &
                      maprec, matrAsse, &
                      faccvg, ldccvg)
        use NonLin_Datastructure_type
        use NonLinearDyna_type
        use Rom_Datastructure_type
        integer(kind=8), intent(in) :: listFuncActi(*)
        character(len=*), intent(in) :: modelz
        character(len=24), intent(in) :: caraElem
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        character(len=19), intent(in) :: listLoad, sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: numeTime
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        type(NL_DS_Contact), intent(inout) :: ds_contact
        type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
        type(NL_DS_Print), intent(inout) :: ds_print
        type(NL_DS_Measure), intent(inout) :: ds_measure
        character(len=19), intent(in) :: hval_algo(*), hval_incr(*)
        character(len=19), intent(in) :: hval_meelem(*), hval_measse(*)
        character(len=24), intent(inout) :: numeDof
        character(len=24), intent(in) :: numeDofFixe
        character(len=19), intent(in) :: solveu
        type(NL_DS_System), intent(in) :: ds_system
        character(len=19), intent(in) :: maprec
        character(len=19), intent(inout) :: matrAsse
        integer(kind=8), intent(out) :: faccvg, ldccvg
    end subroutine nmprma
end interface
