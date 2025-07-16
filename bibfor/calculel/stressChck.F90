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
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.a
! --------------------------------------------------------------------
!
subroutine stressChck(stressModelZ, stressZ, projectOnLigrel, iret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/compareFieldShape.h"
!
    character(len=*), intent(in) :: stressModelZ, stressZ
    aster_logical, intent(in) :: projectOnLigrel
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! Check stress field
!
! --------------------------------------------------------------------------------------------------
!
! In  stressModel      : stress field as model
! In  stress           : stress field
! In  projectOnLigrel  : project field on FED from model
! Out iret             : error code
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
    call compareFieldShape(stressModelZ, stressZ, &
                           projectOnLigrel, "PSIEF_R", &
                           iret)
!
end subroutine
