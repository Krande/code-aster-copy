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
    subroutine nmflma(matrType, mod45,&
                  l_hpp, lModiRigi,&
                  listFuncActi, ds_algopara,&
                  modelZ, caraElem,&
                  ds_material, ds_constitutive,&
                  sddyna, listLoad,&
                  sddisc, numeTime,&
                  ds_posttimestep, nbDofExcl,&
                  hval_incr, hval_algo, &
                  numeDof, ds_system,&
                  ds_measure, hval_meelem,&
                  matrAsse, matrGeom)
        use NonLin_Datastructure_type
        character(len=16), intent(in) :: matrType
        character(len=4), intent(in) :: mod45
        aster_logical, intent(in) :: l_hpp, lModiRigi
        integer(kind=8), intent(in) :: listFuncActi(*)
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        character(len=*), intent(in) :: modelZ
        character(len=24), intent(in) :: caraElem
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        character(len=19), intent(in) :: sddyna, listLoad
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: numeTime
        type(NL_DS_PostTimeStep), intent(in) :: ds_posttimestep
        integer(kind=8), intent(in) :: nbDofExcl
        character(len=19), intent(in) :: hval_algo(*), hval_incr(*)
        character(len=24), intent(in) :: numeDof
        type(NL_DS_System), intent(in) :: ds_system
        type(NL_DS_Measure), intent(inout) :: ds_measure
        character(len=19), intent(in) :: hval_meelem(*)
        character(len=19), intent(out) :: matrAsse, matrGeom
    end subroutine nmflma
end interface
