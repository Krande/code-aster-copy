! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine carc_chck(behaviourPrepCrit)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/utmess.h"
!
    type(Behaviour_PrepCrit), intent(in) :: behaviourPrepCrit
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Some checks
!
! --------------------------------------------------------------------------------------------------
!
! In  behaviourPrepCrit: datastructure to prepare parameters for constitutive laws
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iFactorKeyword, nbFactorKeyword
    aster_logical :: l_mfront_proto, l_mfront_offi, l_umat, lProtoAQ
    character(len=16) :: rela_comp
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = behaviourPrepCrit%nb_comp
    lProtoAQ = ASTER_FALSE

! - Loop on occurrences of COMPORTEMENT
    do iFactorKeyword = 1, nbFactorKeyword
        rela_comp = behaviourPrepCrit%v_crit(iFactorKeyword)%rela_comp

! ----- Detection of specific cases
        call comp_meca_l(rela_comp, 'MFRONT_PROTO', l_mfront_proto)
        call comp_meca_l(rela_comp, 'MFRONT_OFFI', l_mfront_offi)
        call comp_meca_l(rela_comp, 'UMAT', l_umat)

        if (l_mfront_proto .or. l_umat) then
            lProtoAQ = ASTER_TRUE
        end if

! ----- Ban if RELATION = MFRONT and ITER_INTE_PAS negative
        if (behaviourPrepCrit%v_crit(iFactorKeyword)%iter_inte_pas .lt. 0.d0) then
            if (l_mfront_offi .or. l_mfront_proto) then
                call utmess('F', 'COMPOR1_95')
            end if
        end if
    end do

    if (lProtoAQ) then
        call utmess('A', 'QUALITY1_3')
    end if
!
end subroutine
