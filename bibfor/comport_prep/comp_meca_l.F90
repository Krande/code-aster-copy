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
recursive subroutine comp_meca_l(rela_comp, whatz, l_detec, post_iter)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterc/lccree.h"
#include "asterc/lctype.h"
#include "asterc/lcdiscard.h"
!
    character(len=16), intent(in) :: rela_comp
    character(len=*), intent(in) :: whatz
    aster_logical, intent(out) :: l_detec
    character(len=16), optional, intent(in) :: post_iter
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Detection of specific cases
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp    : RELATION comportment
! In  what         : what to detect
! In  post_iter    : type of post_treatment
! Out l_detec      : .true. if specific case is detected
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: what, rela_comp_py, ldctyp
!
! --------------------------------------------------------------------------------------------------
!
    what = whatz
    l_detec = .false.
    if (what .eq. 'CRISTAL') then
        l_detec = (rela_comp .eq. 'MONOCRISTAL') .or. (rela_comp .eq. 'POLYCRISTAL')
    else if (what .eq. 'KIT_META') then
        l_detec = (rela_comp .eq. 'KIT_META')
    else if (what .eq. 'KIT_THM') then
        l_detec = ((rela_comp(1:5) .eq. 'KIT_H') .or. (rela_comp(1:6) .eq. 'KIT_TH'))
    else if (what .eq. 'KIT_DDI') then
        l_detec = (rela_comp .eq. 'KIT_DDI')
    else if (what .eq. 'KIT_CG') then
        l_detec = (rela_comp .eq. 'KIT_CG')
    else if (what .eq. 'KIT') then
        l_detec = (rela_comp(1:4) .eq. 'KIT_')
    else if (what .eq. 'JOINT_MECA_FROT') then
        l_detec = (rela_comp .eq. 'JOINT_MECA_FROT')
    else if (what .eq. 'JOINT_MECA_RUPT') then
        l_detec = (rela_comp .eq. 'JOINT_MECA_RUPT')
    else if (what .eq. 'JOINT_MECA_ENDO') then
        l_detec = (rela_comp .eq. 'JOINT_MECA_ENDO')
    else if (what .eq. 'UMAT') then
        l_detec = (rela_comp .eq. 'UMAT')
    else if (what .eq. 'MFRONT_OFFI') then
        call lccree(1, rela_comp, rela_comp_py)
        call lctype(rela_comp_py, ldctyp)
        call lcdiscard(rela_comp_py)
        l_detec = ldctyp == 'mfront'
    else if (what .eq. 'MFRONT_PROTO') then
        l_detec = (rela_comp .eq. 'MFRONT')
    else if (what .eq. 'MFRONT') then
        call comp_meca_l(rela_comp, 'MFRONT_PROTO', l_detec)
        if (.not. l_detec) then
            call comp_meca_l(rela_comp, 'MFRONT_OFFI', l_detec)
        end if
    else if (what .eq. 'EXTE_COMP') then
        call comp_meca_l(rela_comp, 'MFRONT', l_detec)
        if (.not. l_detec) then
            call comp_meca_l(rela_comp, 'UMAT', l_detec)
        end if
    else if (what .eq. 'PMF') then
        l_detec = (rela_comp .eq. 'MULTIFIBRE')
    else if (what .eq. 'CRIT_RUPT') then
        ASSERT(present(post_iter))
        l_detec = post_iter .eq. 'CRIT_RUPT'
    else
        write (6, *) 'What: ', rela_comp, '/', what, '/', whatz
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
