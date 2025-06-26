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
subroutine comp_meta_info(factorKeyword, metaPrepBehaviour)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/getfac.h"
!
    character(len=16), intent(in) :: factorKeyword
    type(META_PrepBehaviour), intent(out) :: metaPrepBehaviour
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Create datastructure to prepare parameters for behaviour of metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read
! Out metaPrepBehaviour: datastructure to prepare parameters for behaviour of metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_info_comp, nbFactorKeyword
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = 0
    metaPrepBehaviour%factorKeyword = factorKeyword

! - Number of behaviours
    call getfac(factorKeyword, nbFactorKeyword)
    if (nbFactorKeyword .eq. 0) then
        nb_info_comp = 1
    else
        nb_info_comp = nbFactorKeyword
    end if

! - Save number of behaviours
    metaPrepBehaviour%nbFactorKeyword = nbFactorKeyword

! - Allocate comportment informations objects
    allocate (metaPrepBehaviour%paraBehaviour(nb_info_comp))
!
end subroutine
