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
    subroutine nmener(valinc, veasse, measse,&
                      sddyna, nlDynaDamping,&
                      eta, ds_energy, listFuncActi, numeDof, numeDofFixe,&
                      meelem, numeTime, model, ds_material, caraElem,&
                      ds_constitutive, ds_measure, sddisc, solalg,&
                      ds_contact, ds_system)
        use NonLin_Datastructure_type
        use NonLinearDyna_type
        character(len=19), intent(in) :: valinc(*), veasse(*), measse(*)
        character(len=19), intent(in) :: sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        type(NL_DS_Energy), intent(inout) :: ds_energy
        type(NL_DS_Material), intent(in) :: ds_material
        character(len=19), intent(in) :: meelem(*), sddisc, solalg(*)
        character(len=24), intent(in) :: numeDof, numeDofFixe, model, caraElem
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_System), intent(in) :: ds_system
        type(NL_DS_Measure), intent(inout) :: ds_measure
        real(kind=8), intent(in) :: eta
        integer(kind=8), intent(in) :: listFuncActi(*), numeTime
        type(NL_DS_Contact), intent(in) :: ds_contact
    end subroutine nmener
end interface
