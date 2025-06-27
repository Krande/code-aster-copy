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
    subroutine nmflal(optionSpec, ds_posttimestep,&
                  mod45, l_hpp,&
                  nbFreq, coefDimSpace, matrType, optionModal, bande,&
                  nbDofExcl, nbDofStab, lModiRigi, calcLevel)
        use NonLin_Datastructure_type
        character(len=16), intent(in) :: optionSpec
        type(NL_DS_PostTimeStep), intent(in) :: ds_posttimestep
        character(len=16), intent(out) :: optionModal
        character(len=4), intent(out) :: mod45
        integer(kind=8), intent(out) :: nbFreq, nbDofExcl, nbDofStab, coefDimSpace
        character(len=16), intent(out) :: matrType
        aster_logical, intent(out) :: lModiRigi
        real(kind=8), intent(out) :: bande(2)
        aster_logical, intent(out) :: l_hpp
        character(len=16), intent(out) :: calcLevel
    end subroutine nmflal
end interface
