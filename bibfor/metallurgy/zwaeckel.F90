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
! -----------------------------------------2---------------------------
! aslint: disable=W0413
!
subroutine zwaeckel(metaSteelPara, nbPhase, nbVari, &
                    temp0, temp1, temp2, &
                    deltaTime01, deltaTime12, &
                    metaPrev, metaCurr)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/smcarc.h"
#include "asterfort/utmess.h"
#include "asterfort/metaSteelGrainSize.h"
#include "asterfort/Metallurgy_type.h"
!
    type(META_SteelParameters), intent(in) :: metaSteelPara
    integer, intent(in) :: nbPhase, nbVari
    real(kind=8), intent(in) :: temp0, temp1, temp2
    real(kind=8), intent(in) :: deltaTime01, deltaTime12
    real(kind=8), intent(in) :: metaPrev(nbVari)
    real(kind=8), intent(out) :: metaCurr(nbVari)
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY -  Waeckel law for steel phase computing
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
    real(kind=8) :: vari_dumm(9)
    real(kind=8) :: dmoins
    real(kind=8) :: tpoint, ti, tpi
    real(kind=8) :: zeq1, zeq2, zaust, z2, epsi, dt21_mod
    real(kind=8) :: coef_phase
    real(kind=8) :: zeq1i, zeq2i, ti1, ti2, taux
    integer :: i, j, nbpas, nbHist
    aster_logical :: l_cold
    real(kind=8), parameter :: un = 1.d0, zero = 0.d0
!
! --------------------------------------------------------------------------------------------------
!
    epsi = 1.d-10
    ASSERT(nbPhase .eq. NBPHASESTEEL)
    ASSERT(nbVari .eq. NBVARIWAECKEL)

! - Get material parameters for steel
    nbHist = metaSteelPara%trc%nbHist

! - Temperature
    metaCurr(nbPhase+STEEL_TEMP) = temp2
    metaCurr(nbPhase+TEMP_MARTENSITE) = metaPrev(nbPhase+TEMP_MARTENSITE)
    tpoint = (temp1-temp0)/deltaTime01

! - Proportion of austenite
    zaust = un-(metaPrev(PFERRITE)+metaPrev(PPERLITE)+ &
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
    else if (zaust .le. zeq2) then
        l_cold = ASTER_FALSE
    else
        l_cold = ASTER_TRUE
    end if
!
    if (l_cold) then
        if (abs(temp2-temp1) .gt. 5.001d0) then
            nbpas = int(abs(temp2-temp1)/5.d0-0.001d0)+1
            dt21_mod = deltaTime12/dble(nbpas)
            vari_dumm(:) = metaPrev(:)
            do i = 1, nbpas
                ti = temp1+(temp2-temp1)*dble(i-1)/dble(nbpas)
                metaCurr(nbPhase+STEEL_TEMP) = temp1+(dble(i)*(temp2-temp1))/dble(nbpas)
                tpi = (metaCurr(nbPhase+STEEL_TEMP)-ti)/dt21_mod
                call smcarc(nbHist, nbPhase, &
                            zr(metaSteelPara%trc%jv_ftrc), &
                            zr(metaSteelPara%trc%jv_trc), &
                            zr(metaSteelPara%trc%iadtrc+3), &
                            zr(metaSteelPara%trc%iadtrc+metaSteelPara%trc%iadexp), &
                            metaSteelPara, &
                            ti, tpi, deltaTime01, &
                            vari_dumm, metaCurr)
                vari_dumm(:) = metaCurr(:)
            end do
        else
            call smcarc(nbHist, nbPhase, &
                        zr(metaSteelPara%trc%jv_ftrc), &
                        zr(metaSteelPara%trc%jv_trc), &
                        zr(metaSteelPara%trc%iadtrc+3), &
                        zr(metaSteelPara%trc%iadtrc+metaSteelPara%trc%iadexp), &
                        metaSteelPara, &
                        temp1, tpoint, deltaTime01, &
                        metaPrev, metaCurr)
        end if
    else
        if (abs(temp2-temp1) .gt. 5.001d0) then
! ----------------SUBDIVISION EN PAS DE CING DEGRE MAX
            nbpas = int(abs(temp2-temp1)/5.d0-0.001d0)+1
            dt21_mod = deltaTime12/dble(nbpas)
            dmoins = metaPrev(nbPhase+SIZE_GRAIN)
            do i = 1, nbpas
                ti1 = temp1+(temp2-temp1)*dble(i-1)/dble(nbpas)
                ti2 = temp1+(temp2-temp1)*dble(i)/dble(nbpas)
                tpoint = (ti2-ti1)/dt21_mod
                zeq1i = min((ti1-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
                zeq2i = min((ti2-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
                taux = metaSteelPara%taux_1+(metaSteelPara%taux_3-metaSteelPara%taux_1)*zeq1i
                if ((ti1 .lt. (metaSteelPara%ac1-epsi)) .or. (zaust .ge. un)) then
                    z2 = zaust
                else
                    if (zeq2i .ge. (un-epsi)) then
                        tpoint = zero
                    end if
                    if (zaust .gt. zeq1i) then
                        z2 = (taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))
                        z2 = z2*exp(-dt21_mod/taux)
                        z2 = ((-taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))+ &
                              zeq2i+z2-zeq1i)* &
                             (un-zaust)/(un-zeq1i)
                        z2 = z2+zaust
                    else
                        z2 = (taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))-zeq1i+zaust
                        z2 = z2*exp(-dt21_mod/taux)
                        z2 = (-taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))+zeq2i+z2
                    end if
                end if
! ------------- Compute size of grain
                coef_phase = 1.d0
                if (abs(z2) .ge. epsi) then
                    coef_phase = (1.d0-(z2-zaust)/z2)
                end if
                call metaSteelGrainSize(metaSteelPara, &
                                        ti1, deltaTime01, deltaTime12, &
                                        z2, coef_phase, &
                                        dmoins, metaCurr(nbPhase+SIZE_GRAIN))
                if (metaSteelPara%l_grain_size) then
                    zaust = z2
                    dmoins = metaCurr(nbPhase+SIZE_GRAIN)
                end if
            end do
        else
            dt21_mod = deltaTime12
            taux = metaSteelPara%taux_1+(metaSteelPara%taux_3-metaSteelPara%taux_1)*zeq1
            if ((temp1 .lt. (metaSteelPara%ac1-epsi)) .or. (zaust .ge. un)) then
                z2 = zaust
            else
                if (zeq2 .ge. (un-epsi)) then
                    tpoint = zero
                end if
                if (zaust .gt. zeq1) then
                    z2 = (taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))
                    z2 = z2*exp(-dt21_mod/taux)
                    z2 = ((-taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))+zeq2+z2-zeq1)* &
                         (un-zaust)/(un-zeq1)
                    z2 = z2+zaust
                else
                    z2 = (taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))-zeq1+zaust
                    z2 = z2*exp(-dt21_mod/taux)
                    z2 = (-taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))+zeq2+z2
                end if
            end if
! --------- Compute size of grain
            coef_phase = 1.d0
            if (abs(z2) .ge. epsi) then
                coef_phase = (1.d0-(z2-zaust)/z2)
            end if
            call metaSteelGrainSize(metaSteelPara, &
                                    temp1, deltaTime01, deltaTime12, &
                                    z2, coef_phase, &
                                    metaPrev(nbPhase+SIZE_GRAIN), metaCurr(nbPhase+SIZE_GRAIN))
        end if
!           REPARTITION DE DZGAMMA SUR DZALPHA
        if (z2 .gt. (un-epsi)) then
            z2 = un
        end if
        if (zaust .ne. un) then
            do j = 1, 4
                metaCurr(j) = metaPrev(j)*(un-(z2-zaust)/(un-zaust))
            end do
        else
            do j = 1, 4
                metaCurr(j) = metaPrev(j)
            end do
        end if
    end if

! - Compute austenite
    zaust = un-metaCurr(1)-metaCurr(2)-metaCurr(3)-metaCurr(4)
    if (zaust .le. r8prem()) then
        zaust = 0.d0
    end if
    if (zaust .ge. 1.d0) then
        zaust = 1.d0
    end if
    metaCurr(PAUSTENITE) = zaust

! - Compute sum of cold phases
    metaCurr(PSUMCOLD) = 1.d0-zaust
!
end subroutine
