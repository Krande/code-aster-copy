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
    subroutine nmprde(modele, numedd, numfix, ds_material, carele, &
                      ds_constitutive, lischa, ds_algopara, solveu, ds_system, &
                      fonact, ds_print, ds_measure, ds_algorom, sddisc, numins, &
                      valinc, solalg, matass, maprec, ds_contact, &
                      sddyna, nlDynaDamping, &
                      meelem, measse, veelem, veasse, &
                      ldccvg, faccvg, rescvg)
        use NonLin_Datastructure_type
        use Rom_Datastructure_type
        use NonLinearDyna_type
        integer(kind=8) :: fonact(*)
        integer(kind=8) :: numins, ldccvg, faccvg, rescvg
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        character(len=19) :: maprec, matass
        type(NL_DS_Measure), intent(inout) :: ds_measure
        type(NL_DS_Print), intent(inout) :: ds_print
        type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
        character(len=19) :: lischa, solveu, sddisc
        character(len=19), intent(in) :: sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        character(len=24) :: numedd, numfix
        character(len=24) :: modele, carele
        type(NL_DS_System), intent(in) :: ds_system
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_Contact), intent(inout) :: ds_contact
        character(len=19) :: veelem(*), veasse(*)
        character(len=19) :: meelem(*), measse(*)
        character(len=19) :: solalg(*), valinc(*)
    end subroutine nmprde
end interface
