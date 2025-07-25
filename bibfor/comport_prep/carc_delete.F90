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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine carc_delete(prepMapCarcri)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
!
    type(BehaviourPrep_MapCarcri), intent(inout) :: prepMapCarcri
!
! --------------------------------------------------------------------------------------------------
!
! Parameters for integration of constitutive laws (mechanics)
!
! Delete objects
!
! --------------------------------------------------------------------------------------------------
!
! IO  prepMapCarcri    : datastructure to construct CARCRI map
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbFactorKeyword, nbInfo, iInfo
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = prepMapCarcri%nb_comp
    if (nbFactorKeyword .eq. 0) then
        nbInfo = 1
    else
        nbInfo = nbFactorKeyword
    end if
    do iInfo = 1, nbInfo
        if (associated(prepMapCarcri%prepCrit(iInfo)%resi_inte)) then
            deallocate (prepMapCarcri%prepCrit(iInfo)%resi_inte)
        end if
        if (associated(prepMapCarcri%prepCrit(iInfo)%iter_inte_maxi)) then
            deallocate (prepMapCarcri%prepCrit(iInfo)%iter_inte_maxi)
        end if
    end do
    deallocate (prepMapCarcri%prepCrit)
!
end subroutine
