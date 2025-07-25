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
    subroutine nmdesc(modele, numedd, &
                      numfix, ds_material, carele, &
                      ds_constitutive, lischa, ds_contact, &
                      ds_algopara, ds_system, solveu, &
                      fonact, numins, iterat, &
                      sddisc, ds_print, ds_measure, &
                      ds_algorom, sddyna, nlDynaDamping, sdnume, &
                      sderro, matass, maprec, &
                      valinc, solalg, meelem, &
                      measse, veasse, lerrit)
        use NonLin_Datastructure_type
        use NonLinearDyna_type
        use ROM_Datastructure_type
        use HHO_type
        character(len=24) :: modele
        character(len=24) :: numedd
        character(len=24) :: numfix
        type(NL_DS_Material), intent(in) :: ds_material
        character(len=24) :: carele
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        character(len=19), intent(in) :: lischa, sddyna
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        type(NL_DS_Contact), intent(inout) :: ds_contact
        character(len=19) :: solveu
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        integer(kind=8) :: fonact(*)
        integer(kind=8) :: numins
        integer(kind=8) :: iterat
        character(len=19) :: sddisc
        type(NL_DS_Print), intent(inout) :: ds_print
        type(NL_DS_Measure), intent(inout) :: ds_measure
        type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
        type(NL_DS_System), intent(in) :: ds_system
        character(len=19) :: sdnume
        character(len=24) :: sderro
        character(len=19) :: matass
        character(len=19) :: maprec
        character(len=19) :: valinc(*)
        character(len=19) :: solalg(*)
        character(len=19) :: meelem(*)
        character(len=19) :: measse(*)
        character(len=19) :: veasse(*)
        aster_logical :: lerrit
    end subroutine nmdesc
end interface
