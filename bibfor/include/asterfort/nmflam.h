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
    subroutine nmflam(optionSpec,&
                  model, ds_material, caraElem, listLoad, listFuncActi,&
                  numeDof, ds_system ,&
                  ds_constitutive, &
                  sddisc, numeTime,&
                  sddyna, sderro, ds_algopara,&
                  ds_measure,&
                  hval_incr, hval_algo,&
                  hval_meelem, &
                  ds_posttimestep)
        use NonLin_Datastructure_type
        character(len=16), intent(in) :: optionSpec
        character(len=24), intent(in) :: model, caraElem
        type(NL_DS_Material), intent(in) :: ds_material
        character(len=19), intent(in) :: listLoad
        integer(kind=8), intent(in) :: listFuncActi(*)
        character(len=24), intent(in) :: numeDof
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: numeTime
        character(len=19), intent(in) :: sddyna
        character(len=24), intent(in) :: sderro
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        type(NL_DS_Measure), intent(inout) :: ds_measure
        type(NL_DS_System), intent(in) :: ds_system
        character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
        character(len=19), intent(in) :: hval_meelem(*)
        type(NL_DS_PostTimeStep), intent(inout) :: ds_posttimestep
    end subroutine nmflam
end interface
