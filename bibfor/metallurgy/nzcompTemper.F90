! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine nzcompTemper(metaParaTemper, numeCompTemper, &
                        nbVari, nbVariTemper, &
                        deltaTime12, &
                        infoTemper, metaIn, metaOut)
!
    use Metallurgy_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "asterfort/zjma.h"
#include "asterfort/Metallurgy_type.h"
!
    type(META_MaterialParameters), intent(in) :: metaParaTemper
    integer(kind=8), intent(in) :: numeCompTemper
    integer(kind=8), intent(in) :: nbVari, nbVariTemper
    real(kind=8), intent(in) :: deltaTime12
    real(kind=8), intent(in) :: infoTemper(NB_PARAIN_TEMPER)
    real(kind=8), intent(in) :: metaIn(nbVari)
    real(kind=8), intent(out) :: metaOut(nbVariTemper)
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY - Compute phases (tempering case)
!
! General
!
! --------------------------------------------------------------------------------------------------
!
! In  metaParaTemper      : material parameters for tempering for metallurgy
! In  numeCompTemper      : index of tempering law for metallurgy
! In  nbVari              : number of internal state variables without tempering
! In  nbVariTemper        : number of internal state variables with tempering
! In  deltaTime12         : increment of time [N, N+1]
! In  temp1               : temperature at time N
! In  temp2               : temperature at time N+1
! In  infoTemper          : value parameters for tempering
! In  metaIn              : value of internal state variable without tempering
! Out metaOut             : value of internal state variable with tempering
!
! --------------------------------------------------------------------------------------------------
!
    select case (numeCompTemper)

    case (3)
        call zjma(metaParaTemper%steel, &
                  nbVari, nbVariTemper, &
                  deltaTime12, &
                  infoTemper, metaIn, metaOut)

    case default
        call utmess('F', 'COMPOR1_43', si=numeCompTemper)

    end select
!
end subroutine
