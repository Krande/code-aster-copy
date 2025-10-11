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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine comp_ntvari(ligrel_, comporMap_, comporList_, comporInfo, &
                           nt_vari, nb_vari_maxi, mapNbZone, prepExte)
        use BehaviourPrepare_type
        character(len=19), optional, intent(in) :: ligrel_
        character(len=19), optional, intent(in) :: comporMap_
        character(len=16), optional, intent(in) :: comporList_(COMPOR_SIZE)
        character(len=19), intent(in) :: comporInfo
        integer(kind=8), intent(out) :: nt_vari, nb_vari_maxi, mapNbZone
        type(BehaviourPrep_Exte), pointer :: prepExte(:)
    end subroutine comp_ntvari
end interface
