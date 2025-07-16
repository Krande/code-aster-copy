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
subroutine zacier(metaSteelPara, nbPhase, nbVari, &
                  temp0, temp1, temp2, &
                  deltaTime01, deltaTime12, &
                  metaPrev, metaCurr)
!
    use Metallurgy_type
    use MetallurgySteel_Compute_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/zwaeckel.h"
!
    type(META_SteelParameters), intent(in) :: metaSteelPara
    integer(kind=8), intent(in) :: nbPhase, nbVari
    real(kind=8), intent(in) :: temp0, temp1, temp2
    real(kind=8), intent(in) :: deltaTime01, deltaTime12
    real(kind=8), intent(in) :: metaPrev(nbVari)
    real(kind=8), intent(out) :: metaCurr(nbVari)
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY -  Compute phase
!
! Main law for steel
!
! --------------------------------------------------------------------------------------------------
!
! In  metaSteelPara       : material parameters for metallurgy of steel
! In  nbPhase             : number of phases
! In  nbVari              : number of internal state variables
! In  temp0               : temperature at time N-1
! In  temp1               : temperature at time N
! In  temp2               : temperature at time N+1
! In  deltaTime01         : increment of time [N-1, N]
! In  deltaTime12         : increment of time [N, N+1]
! In  metaPrev            : value of internal state variable at previous time step
! Out metaCurr            : value of internal state variable at current time step
!
! --------------------------------------------------------------------------------------------------
!
    call zwaeckel(metaSteelPara, nbPhase, nbVari, &
                  temp0, temp1, temp2, &
                  deltaTime01, deltaTime12, &
                  metaPrev, metaCurr)
    call metaSteelCheckPhases(nbVari, metaCurr)
!
end subroutine
