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
subroutine nzcomp(jvMaterCode, metaPara, numeComp, &
                  nbPhase, nbVari, &
                  deltaTime01, deltaTime12, time2, &
                  tempInit, temp1, temp2, &
                  metaPrev, metaCurr, &
                  lNodeDebug_)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "asterfort/zacier.h"
#include "asterfort/zedgar.h"
#include "asterfort/Metallurgy_type.h"
!
    integer(kind=8), intent(in) :: jvMaterCode
    type(META_MaterialParameters), intent(inout) :: metaPara
    integer(kind=8), intent(in) :: numeComp, nbPhase, nbVari
    real(kind=8), intent(in) :: deltaTime01, deltaTime12, time2
    real(kind=8), intent(in) :: tempInit, temp1, temp2
    real(kind=8), intent(in) :: metaPrev(nbVari)
    real(kind=8), intent(out) :: metaCurr(nbVari)
    aster_logical, optional, intent(in) :: lNodeDebug_
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY - Compute phases
!
! General
!
! --------------------------------------------------------------------------------------------------
!
! In  jvMaterCode         : coded material address
! In  metaPara            : material parameters for metallurgy
! In  numeComp            : index of behaviour law for metallurgy
! In  nbPhase             : number of phases
! In  nbVari              : number of internal state variables
! In  tempInit            : temperature at time N-1
! In  temp1               : temperature at time N
! In  temp2               : temperature at time N+1
! In  deltaTime01         : increment of time [N-1, N]
! In  deltaTime12         : increment of time [N, N+1]
! In  time2               : value of time N+1
! In  metaPrev            : value of internal state variable at previous time step
! Out metaCurr            : value of internal state variable at current time step
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lNodeDebug
!
! --------------------------------------------------------------------------------------------------
!
    lNodeDebug = ASTER_FALSE
    if (present(lNodeDebug_)) then
        lNodeDebug = lNodeDebug_
    end if

    select case (numeComp)
!
    case (0)
! ----- Empty behaviour
    case (1)
        call zedgar(jvMaterCode, nbPhase, &
                    temp1, temp2, &
                    time2, deltaTime12, &
                    metaPrev, metaCurr)
    case (2)
        metaPara%steel%lNodeDebug = lNodeDebug
        call zacier(metaPara%steel, nbPhase, nbVari, &
                    tempInit, temp1, temp2, &
                    deltaTime01, deltaTime12, &
                    metaPrev, metaCurr)
    case default
        call utmess('F', 'COMPOR1_43', si=numeComp)

    end select
!
end subroutine
