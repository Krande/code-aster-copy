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
! aslint: disable=W0413
!
subroutine zacier(metaSteelPara, nbPhase, &
                  temp0, temp1, temp2, &
                  deltaTime01, deltaTime12, &
                  metaPrev, metaCurr)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/metaSteelCheckPhases.h"
#include "asterfort/metaSteelGrainSize.h"
#include "asterfort/smcarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(META_SteelParameters), intent(in) :: metaSteelPara
    integer, intent(in) :: nbPhase
    real(kind=8), intent(in) :: temp0, temp1, temp2
    real(kind=8), intent(in) :: deltaTime01, deltaTime12
    real(kind=8), intent(in) :: metaPrev(nbPhase+3)
    real(kind=8), intent(out) :: metaCurr(nbPhase+3)
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
! In  temp0               : temperature at time N-1
! In  temp1               : temperature at time N
! In  temp2               : temperature at time N+1
! In  deltaTime01         : increment of time [N-1, N]
! In  deltaTime12         : increment of time [N, N+1]
! In  metaPrev            : value of internal state variable at previous time step
! Out metaPrev            : value of internal state variable at current time step
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: epsi = 1.d-10, epsiTemp = 0.001d0
    real(kind=8) :: metaDumm(nbPhase+3)
    real(kind=8) :: dPrev, dStep, dCurr
    real(kind=8) :: tempStep, tPoint
    real(kind=8) :: zeq1, zeq2, zAustPrev, zAustCurr, deltaTimeStep
    real(kind=8) :: coef_phase
    real(kind=8) :: ti1, ti2, taux
    integer :: iStep, j, nbStep, nbHist, nbVari
    aster_logical :: l_cold, lMultiStep
    real(kind=8), parameter :: un = 1.d0, zero = 0.d0
!
! --------------------------------------------------------------------------------------------------
!
    nbVari = nbPhase+3
    ASSERT(nbPhase .eq. 6)

! - Get material parameters for steel
    nbHist = metaSteelPara%trc%nbHist

! - Temperature
    metaCurr(nbPhase+STEEL_TEMP) = temp2
    metaCurr(nbPhase+TEMP_MARTENSITE) = metaPrev(nbPhase+TEMP_MARTENSITE)
    tpoint = (temp1-temp0)/deltaTime01

! - Proportion of austenite
    zAustPrev = un-(metaPrev(PFERRITE)+metaPrev(PPERLITE)+ &
                    metaPrev(PBAINITE)+metaPrev(PMARTENS))

! - Colding or not ?
    zeq1 = min((temp1-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
    zeq2 = min((temp2-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
    if (tpoint .gt. zero) then
        l_cold = ASTER_FALSE
    else if (temp2 .gt. metaSteelPara%ar3) then
        l_cold = ASTER_FALSE
    else if (temp2 .lt. metaSteelPara%ac1) then
        l_cold = ASTER_TRUE
    else if (tpoint .lt. zero) then
        l_cold = ASTER_TRUE
    else if (zAustPrev .le. zeq2) then
        l_cold = ASTER_FALSE
    else
        l_cold = ASTER_TRUE
    end if

! - Initial number of steps
    lMultiStep = ASTER_FALSE
    nbStep = 1
    if (abs(temp2-temp1) .gt. STEEL_MAX_TEMP_STEP+epsiTemp) then
        nbStep = int(abs(temp2-temp1)/5.d0-epsiTemp)+1
    end if
!
    if (l_cold) then
        deltaTimeStep = deltaTime12/dble(nbStep)
        metaDumm = metaPrev
        do iStep = 1, nbStep
            tempStep = temp1+(temp2-temp1)*dble(iStep-1)/dble(nbStep)
            metaCurr(nbPhase+STEEL_TEMP) = temp1+ &
                                           (dble(iStep)*(temp2-temp1))/dble(nbStep)
            if (nbStep .eq. 1) then
                tPoint = (temp1-temp0)/deltaTime01
            else
                tPoint = (metaCurr(nbPhase+STEEL_TEMP)-tempStep)/deltaTimeStep
            end if
            call smcarc(nbHist, nbPhase, &
                        zr(metaSteelPara%trc%jv_ftrc), &
                        zr(metaSteelPara%trc%jv_trc), &
                        zr(metaSteelPara%trc%iadtrc+3), &
                        zr(metaSteelPara%trc%iadtrc+metaSteelPara%trc%iadexp), &
                        metaSteelPara, &
                        tempStep, tPoint, deltaTime01, &
                        metaDumm, metaCurr)
            metaDumm(:) = metaCurr(:)
        end do
    else
10      continue

        dPrev = metaPrev(nbPhase+SIZE_GRAIN)
        deltaTimeStep = deltaTime12/dble(nbStep)
        dStep = dPrev
        do iStep = 1, nbStep
            if (nbStep .eq. 1) then
                ti1 = temp1
                ti2 = temp2
                tPoint = (temp1-temp0)/deltaTime01
            else
                ti1 = temp1+(temp2-temp1)*dble(iStep-1)/dble(nbStep)
                ti2 = temp1+(temp2-temp1)*dble(iStep)/dble(nbStep)
                tpoint = (ti2-ti1)/deltaTimeStep
            end if
            zeq1 = min((ti1-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
            zeq2 = min((ti2-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
            taux = metaSteelPara%taux_1+(metaSteelPara%taux_3-metaSteelPara%taux_1)*zeq1

            if ((ti1 .lt. (metaSteelPara%ac1-epsi)) .or. (zAustPrev .ge. un)) then
                zAustCurr = zAustPrev
            else
                if (zeq2 .ge. (un-epsi)) then
                    tpoint = zero
                end if
                if (zAustPrev .gt. zeq1) then
                    zAustCurr = (taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))
                    zAustCurr = zAustCurr*exp(-deltaTimeStep/taux)
                    zAustCurr = ((-taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))+ &
                                 zeq2+zAustCurr-zeq1)* &
                                (un-zAustPrev)/(un-zeq1)
                    zAustCurr = zAustCurr+zAustPrev
                else
                    zAustCurr = (taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))-zeq1+zAustPrev
                    zAustCurr = zAustCurr*exp(-deltaTimeStep/taux)
                    zAustCurr = (-taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))+zeq2+zAustCurr
                end if
            end if

! --------- Compute size of grain
            coef_phase = 1.d0
            if (abs(zAustCurr) .ge. epsi) then
                coef_phase = (1.d0-(zAustCurr-zAustPrev)/zAustCurr)
            end if
            call metaSteelGrainSize(metaSteelPara, &
                                    ti1, deltaTime01, deltaTime12, &
                                    zAustCurr, coef_phase, &
                                    dPrev, dCurr)
            metaCurr(nbPhase+SIZE_GRAIN) = dCurr
            if (metaSteelPara%l_grain_size .and. nbStep .gt. 1) then
                zAustPrev = zAustCurr
                dStep = dCurr
            end if
        end do

! ----- Save size of grain
        metaCurr(nbPhase+SIZE_GRAIN) = dCurr

! ----- Check bounds
        if (zAustCurr .gt. (un-epsi)) then
            zAustCurr = un
        end if
        if (zAustCurr .lt. 0.d0 .and. zAustCurr .ge. -STEEL_MIN_CUT) then
            zAustCurr = 0.d0
        end if
        if (zAustCurr .lt. -STEEL_MIN_CUT) then
            nbStep = STEEL_MAX_NB_STEP
            if (lMultiStep) then
                call utmess('F', 'META1_51')
            else
                lMultiStep = ASTER_TRUE
                goto 10
            end if
        end if

! ----- Compute new ferrite phases
        if (zAustPrev .ne. un) then
            do j = 1, 4
                metaCurr(j) = metaPrev(j)*(un-(zAustCurr-zAustPrev)/(un-zAustPrev))
            end do
        else
            do j = 1, 4
                metaCurr(j) = metaPrev(j)
            end do
        end if
    end if

! - Compute austenite
    zAustCurr = un-metaCurr(1)-metaCurr(2)-metaCurr(3)-metaCurr(4)
    if (zAustCurr .le. r8prem()) then
        zAustCurr = 0.d0
    end if
    if (zAustCurr .ge. 1.d0) then
        zAustCurr = 1.d0
    end if
    metaCurr(PAUSTENITE) = zAustCurr

! - Compute sum of cold phases
    metaCurr(PSUMCOLD) = 1.d0-zAustCurr

! - Check
    call metaSteelCheckPhases(nbVari, metaCurr)
!
end subroutine
