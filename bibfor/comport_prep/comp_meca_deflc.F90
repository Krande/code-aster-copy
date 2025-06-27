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
! person_in_charge: astrid.filiot at edf.fr
!
subroutine comp_meca_deflc(rela_comp, defo_comp, defo_ldc)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/lccree.h"
#include "asterc/lcdeformldc.h"
#include "asterc/lcdiscard.h"
#include "asterfort/assert.h"
!
    character(len=16), intent(in) :: rela_comp
    character(len=16), intent(in) :: defo_comp
    character(len=16), intent(out) :: defo_ldc
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Select type of strain (mechanical, total or old)
!
! --------------------------------------------------------------------------------------------------
!
! In  rela_comp   : comportement RELATION
! In  defo_comp   : type of deformation (kinematik)
! Out defo_ldc    : type of strain ('MECANIQUE','TOTALE' or 'OLD')
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: rela_code_py, defo_code_py
    character(len=16) :: defo_ldc_rela, defo_ldc_defo
    aster_logical :: l_deformlc
!
! --------------------------------------------------------------------------------------------------
!
! Attribute "deform_ldc" from behaviour catalog can be 'MECANIQUE', 'TOTALE' or 'OLD' only
!
    call lccree(1, rela_comp, rela_code_py)
    call lcdeformldc(rela_code_py, defo_ldc_rela)
    call lccree(1, defo_comp, defo_code_py)
    call lcdeformldc(defo_code_py, defo_ldc_defo)

    if (defo_ldc_defo .eq. 'TOTALE') then
        defo_ldc = 'TOTALE'
    else if (defo_ldc_defo .eq. 'MECANIQUE') then
        defo_ldc = defo_ldc_rela
    else
        ASSERT(ASTER_FALSE)
    end if
!
! deform_ldc should not be 'MECANIQUE' if the kinematik is 'SIMO_MIEHE'
! MFRONT is the only exception
!
    l_deformlc = (defo_ldc_rela .eq. 'MECANIQUE') .and. (defo_comp .eq. 'SIMO_MIEHE') &
                 .and. (rela_comp .ne. 'MFRONT')
    ASSERT(.not. l_deformlc)

    ASSERT((defo_ldc .eq. 'MECANIQUE') .or. (defo_ldc .eq. 'TOTALE') .or. (defo_ldc .eq. 'OLD'))
!
    call lcdiscard(rela_code_py)
    call lcdiscard(defo_code_py)
!
end subroutine
