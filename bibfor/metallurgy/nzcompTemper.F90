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
subroutine nzcompTemper(metaPara, numeComp, &
                        nbVari, nbVariTemper, nbVariPrev, &
                        deltaTime12, &
                        temp1, temp2, &
                        prevMetaIsTemper, &
                        metaPrev, metaCurr, metaCurrTemper)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "asterfort/zjma.h"
#include "asterfort/Metallurgy_type.h"
!
    type(META_MaterialParameters), intent(in) :: metaPara
    integer(kind=8), intent(in) :: numeComp, nbVari, nbVariTemper, nbVariPrev
    real(kind=8), intent(in) :: deltaTime12
    real(kind=8), intent(in) :: temp1, temp2
    aster_logical, intent(in) :: prevMetaIsTemper
    real(kind=8), intent(in) :: metaPrev(nbVariPrev)
    real(kind=8), intent(in) :: metaCurr(nbVari)
    real(kind=8), intent(out) :: metaCurrTemper(nbVariTemper)
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY - Compute phases (tempering case)
!
! General
!
! --------------------------------------------------------------------------------------------------
!
! In  metaPara            : material parameters for metallurgy
! In  numeComp            : index of behaviour law for metallurgy
! In  nbPhase             : number of phases
! In  nbVari              : number of internal state variables
! In  tempInit            : temperature at time N-1
! In  temp1               : temperature at time N
! In  temp2               : temperature at time N+1
! In  deltaTime01         : increment of time [N-1, N]
! In  deltaTime12         : increment of time [N, N+1]
! In  metaPrev            : value of internal state variable at previous time step
! In  metaCurr            : value of internal state variable at current time step without tempering
! Out metaCurrTemper      : value of internal state variable at current time step with tempering
!
! --------------------------------------------------------------------------------------------------
!
    select case (numeComp)

    case (3)
        call zjma(metaPara%steel, &
                  nbVari, nbVariTemper, nbVariPrev, &
                  temp1, temp2, &
                  deltaTime12, &
                  prevMetaIsTemper, &
                  metaPrev, metaCurr, metaCurrTemper)

    case default
        call utmess('F', 'COMPOR1_43', si=numeComp)

    end select
!
end subroutine
