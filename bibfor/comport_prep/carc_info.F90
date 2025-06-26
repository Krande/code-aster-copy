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
subroutine carc_info(prepMapCarcri)
!
    use BehaviourPrepare_type
!
    implicit none
!
#include "asterc/getfac.h"
!
    type(BehaviourPrep_MapCarcri), intent(out) :: prepMapCarcri
!
! --------------------------------------------------------------------------------------------------
!
! Parameters for integration of constitutive laws (mechanics)
!
! Create objects
!
! --------------------------------------------------------------------------------------------------
!
! Out prepMapCarcri    : datastructure to construct CARCRI map
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    integer(kind=8) :: nbInfo, nbFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!

! - Number of factor keywords
    nbFactorKeyword = 0
    call getfac(factorKeyword, nbFactorKeyword)

! - Number of factor keyword information
    if (nbFactorKeyword .eq. 0) then
        nbInfo = 1
    else
        nbInfo = nbFactorKeyword
    end if

! - Initializations
    prepMapCarcri%nb_comp = nbFactorKeyword
    prepMapCarcri%prepCrit => null()
    allocate (prepMapCarcri%prepCrit(nbInfo))
!
end subroutine
