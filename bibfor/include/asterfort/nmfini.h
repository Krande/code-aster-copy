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
    subroutine nmfini(sddyna, nlDynaDamping,&
                        valinc, measse, model, ds_material,&
                        caraElem, ds_constitutive, ds_system,&
                        ds_measure, sddisc, numeTime,&
                        solalg, numeDof, listFuncActi)
        use NonLin_Datastructure_type
        use NonLinearDyna_type
        character(len=19), intent(in) :: valinc(*), measse(*)
        character(len=19), intent(in) :: sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        integer(kind=8), intent(in) :: listFuncActi(*)
        character(len=24), intent(in) :: model, caraElem
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_System), intent(in) :: ds_system
        type(NL_DS_Measure), intent(inout) :: ds_measure
        character(len=24), intent(in) :: numeDof
        character(len=19), intent(in) :: sddisc, solalg(*)
        integer(kind=8), intent(in) :: numeTime
    end subroutine nmfini
end interface
