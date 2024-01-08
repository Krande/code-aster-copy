! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine comp_meta_info(metaPrepPara)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
!
    type(META_PrepPara), intent(out) :: metaPrepPara
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Create datastructure to prepare parameters for behaviour of metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! Out metaPrepPara     : datastructure to prepare parameters for behaviour of metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    integer :: nb_info_comp, nbFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = 0
    call getfac(factorKeyword, nbFactorKeyword)
    metaPrepPara%para => null()

! - Number of behaviours
    if (nbFactorKeyword .eq. 0) then
        nb_info_comp = 1
    else
        nb_info_comp = nbFactorKeyword
    end if

! - Save number of behaviours
    metaPrepPara%nb_comp = nbFactorKeyword

! - Allocate comportment informations objects
    allocate (metaPrepPara%para(nb_info_comp))
!
end subroutine
