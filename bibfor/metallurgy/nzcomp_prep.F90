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
subroutine nzcomp_prep(jvMaterCode, metaType, metaPara)
!
    use Metallurgy_type
    use MetallurgySteel_Compute_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
!
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=16), intent(in) :: metaType
    type(META_MaterialParameters), intent(out) :: metaPara
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY - Compute phases
!
! General (preparation)
!
! --------------------------------------------------------------------------------------------------
!
! In  jvMaterCode      : coded material address
! In  metaType         : type of phase
! Out metaPara         : material parameters for metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    type(META_SteelParameters) :: metaSteelPara
!
! --------------------------------------------------------------------------------------------------
!
    if (metaType .eq. 'ACIER') then
! ----- Get material parameters for steel
        call metaSteelGetParameters(jvMaterCode, metaSteelPara)
! ----- Get material parameters for TRC curve
        call metaSteelTRCGetParameters(jvMaterCode, metaSteelPara)
    else if (metaType .eq. 'ACIER_REVENU') then
! ----- Get material parameters for steel
        call metaSteelGetParameters(jvMaterCode, metaSteelPara)
        call metaSteelTemperGetParameters(jvMaterCode, metaSteelPara)
    elseif (metaType .eq. 'ZIRC') then
! ----- Depending on temperature: too early here !
    elseif (metaType .eq. 'VIDE') then
! ----- Nothing
    else
        ASSERT(ASTER_FALSE)
    end if
!
    metaPara%steel = metaSteelPara
!
end subroutine
