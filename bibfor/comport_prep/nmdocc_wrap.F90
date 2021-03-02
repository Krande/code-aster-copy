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
subroutine nmdocc_wrap(model, chmate, compor, etat_init, implex, verbose)
!
use Behaviour_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/nmdocc.h"
#include "asterfort/comp_info.h"
!
character(len=8), intent(in) :: model, chmate
character(len=19), intent(in) :: compor
integer, intent(in) :: etat_init, implex, verbose
!
! --------------------------------------------------------------------------------------------------
!
! Wrapper to nmdocc - Required for aster_logical
!
! Get parameters from COMPORTEMENT keyword and prepare COMPOR <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  chmate           : name of material field
! In  compor           : name of <CARTE> COMPOR
! In  etat_init        : 1 if initial state is defined
! In  implex           : 1 if IMPLEX method
! In  verbose          : 1 for verbose mode
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_etat_init, l_implex, lVerbose
!
! --------------------------------------------------------------------------------------------------
!
    l_etat_init = etat_init .gt. 0
    l_implex    = implex .gt. 0
    lVerbose    = verbose .gt. 0
    call nmdocc(model, chmate, l_etat_init, l_implex, compor)
    if (lVerbose) then
        call comp_info(model, compor)
    endif
!
end subroutine
