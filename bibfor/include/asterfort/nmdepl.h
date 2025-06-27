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
#include "asterf_types.h"
!
interface
    subroutine nmdepl(modele, numedd , ds_material, carele,&
                      ds_constitutive, lischa, fonact, ds_measure, ds_algopara,&
                      noma, numins , iterat, solveu, matass,&
                      sddyna, nlDynaDamping,&
                      sddisc, sdnume, sdpilo, sderro,&
                      ds_contact, valinc, solalg, veelem, veasse,&
                      eta, ds_conv, ds_system, lerrit)
        use NonLin_Datastructure_type
        use NonLinearDyna_type
        integer(kind=8) :: fonact(*)
        integer(kind=8) :: iterat, numins
        real(kind=8) :: eta
        character(len=8) :: noma
        type(NL_DS_Conv), intent(inout) :: ds_conv
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        character(len=19) :: sddisc, sdnume, sdpilo
        character(len=19), intent(in) :: sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        type(NL_DS_Measure), intent(inout) :: ds_measure
        character(len=19) :: lischa, matass, solveu
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        type(NL_DS_Material), intent(in) :: ds_material
        character(len=24) :: modele, numedd, carele
        character(len=24) :: sderro
        character(len=19) :: veelem(*), veasse(*)
        character(len=19) :: solalg(*), valinc(*)
        type(NL_DS_System), intent(in) :: ds_system
        type(NL_DS_Contact), intent(inout) :: ds_contact
        aster_logical :: lerrit
    end subroutine nmdepl
end interface
