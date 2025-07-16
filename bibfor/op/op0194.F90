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
subroutine op0194()
!
    use Metallurgy_type
    use MetallurgyOperator_module
!
    implicit none
!
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
!
! --------------------------------------------------------------------------------------------------
!
! Command: CALC_META
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iOption
    character(len=16) :: option
    type(META_ParaOperator) :: metaParaOperator
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Get parameters of command
    call metaGetParameters(metaParaOperator)

! - Compute options
    do iOption = 1, metaParaOperator%nbOption
        option = metaParaOperator%listOption(iOption)
        if (option .eq. 'META_ELNO') then
            call metaCompPhases(metaParaOperator)
        else
            call metaCompOtherOptions(metaParaOperator, option)
        end if
    end do

! - Cleaning
    call metaDelParameters(metaParaOperator)
!
    call jedema()
end subroutine
