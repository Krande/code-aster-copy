! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine nzcomp(jvMaterCode, metaPara, &
                  numeComp, nbPhase, nbVari, &
                  dt10, dt21, inst2, &
                  tno0, tno1, tno2, &
                  metaPrev, metaCurr)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/zacier.h"
#include "asterfort/zedgar.h"
#include "asterfort/Metallurgy_type.h"
!
    integer, intent(in) :: jvMaterCode
    type(META_MaterialParameters), intent(in) :: metaPara
    integer, intent(in) :: numeComp, nbPhase, nbVari
    real(kind=8), intent(in) :: dt10, dt21, inst2
    real(kind=8), intent(in) :: tno0, tno1, tno2
    real(kind=8), intent(in) :: metaPrev(nbVari)
    real(kind=8), intent(out) :: metaCurr(nbVari)
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
! In  tno0                : temperature at time N-1
! In  tno1                : temperature at time N
! In  tno2                : temperature at time N+1
! In  dt10                : increment of time [N-1, N]
! In  dt21                : increment of time [N, N+1]
! In  inst2               : value of time N+1
! In  metaPrev            : value of internal state variable at previous time step
! Out metaCurr            : value of internal state variable at current time step
!
! --------------------------------------------------------------------------------------------------
!
    select case (numeComp)
!
    case (1)
        call zedgar(jvMaterCode, nbPhase, &
                    tno1, tno2, &
                    inst2, dt21, &
                    metaPrev, metaCurr)
    case (2)
        call zacier(metaPara%steel, nbPhase, nbVari, &
                    tno0, tno1, tno2, &
                    dt10, dt21, &
                    metaPrev, metaCurr)
    case default
        ASSERT(ASTER_FALSE)
    end select
!
end subroutine
