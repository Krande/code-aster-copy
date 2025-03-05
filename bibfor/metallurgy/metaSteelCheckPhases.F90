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
subroutine metaSteelCheckPhases(nbVari, metaCurr)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
    integer, intent(in) :: nbVari
    real(kind=8), intent(in) :: metaCurr(nbVari)
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY
!
! Check domain of definition of phases
!
! --------------------------------------------------------------------------------------------------
!
! In  metaCurr            : value of internal state variable at current time step
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iPhase
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbVari .ge. PSTEEL_NB)
    do iPhase = 1, PSTEEL_NB
        if (metaCurr(iPhase) .lt. 0.d0 .or. metaCurr(iPhase) .gt. 1.d0) then
            call utmess('F', 'META1_52', si=iPhase, sr=metaCurr(iPhase))
        end if
    end do
!
end subroutine
