! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine nmdocr_wrap(model, carcri, implex)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/nmdocr.h"
!
character(len=8), intent(in)   :: model
character(len=24), intent(out) :: carcri
integer, intent(in) :: implex
!
! --------------------------------------------------------------------------------------------------
!
! Wrapper to nmdocr - Required for aster_logical
!
! Get parameters from COMPORTEMENT keyword and prepare CARCRI <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! Out carcri           : name of <CARTE> CARCRI
! In  implex           : 1 if IMPLEX method
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_implex
!
! --------------------------------------------------------------------------------------------------
!
    l_implex = implex .gt. 0
    call nmdocr(model, carcri, l_implex)
!
end subroutine
