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
! aslint: disable=W0413
!
subroutine zwaeckel(metaSteelPara, nbPhase, nbVari, &
                    tpg0, tpg1, tpg2, &
                    dt10, dt21, &
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
    real(kind=8), intent(in) :: tpg0, tpg1, tpg2
    real(kind=8), intent(in) :: dt10, dt21
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
! In  tpg0                : temperature at time N-1
! In  tpg1                : temperature at time N
! In  tpg2                : temperature at time N+1
! In  dt10                : increment of time [N-1, N]
! In  dt21                : increment of time [N, N+1]
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
    ASSERT(nbPhase .eq. 6)
    ASSERT(nbVari .eq. 9)

! - Get material parameters for steel
    nbHist = metaSteelPara%trc%nbHist

! - Temperature
    metaCurr(nbPhase+STEEL_TEMP) = tpg2
    metaCurr(nbPhase+TEMP_MARTENSITE) = metaPrev(nbPhase+TEMP_MARTENSITE)
    tpoint = (tpg1-tpg0)/dt10

! - Proportion of austenite
    zaust = un-(metaPrev(PFERRITE)+metaPrev(PPERLITE)+ &
                metaPrev(PBAINITE)+metaPrev(PMARTENS))

! - Colding or not ?
    zeq1 = min((tpg1-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
    zeq2 = min((tpg2-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
    if (tpoint .gt. zero) then
        l_cold = ASTER_FALSE
    else if (tpg2 .gt. metaSteelPara%ar3) then
        l_cold = ASTER_FALSE
    else if (tpg2 .lt. metaSteelPara%ac1) then
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
        if (abs(tpg2-tpg1) .gt. 5.001d0) then
            nbpas = int(abs(tpg2-tpg1)/5.d0-0.001d0)+1
            dt21_mod = dt21/dble(nbpas)
            vari_dumm(:) = metaPrev(:)
            do i = 1, nbpas
                ti = tpg1+(tpg2-tpg1)*dble(i-1)/dble(nbpas)
                metaCurr(nbPhase+STEEL_TEMP) = tpg1+(dble(i)*(tpg2-tpg1))/dble(nbpas)
                tpi = (metaCurr(nbPhase+STEEL_TEMP)-ti)/dt21_mod
                call smcarc(nbHist, nbPhase, &
                            zr(metaSteelPara%trc%jv_ftrc), &
                            zr(metaSteelPara%trc%jv_trc), &
                            zr(metaSteelPara%trc%iadtrc+3), &
                            zr(metaSteelPara%trc%iadtrc+metaSteelPara%trc%iadexp), &
                            metaSteelPara, &
                            ti, tpi, dt10, &
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
                        tpg1, tpoint, dt10, &
                        metaPrev, metaCurr)
        end if
    else
        if (abs(tpg2-tpg1) .gt. 5.001d0) then
! ----------------SUBDIVISION EN PAS DE CING DEGRE MAX
            nbpas = int(abs(tpg2-tpg1)/5.d0-0.001d0)+1
            dt21_mod = dt21/dble(nbpas)
            dmoins = metaPrev(nbPhase+SIZE_GRAIN)
            do i = 1, nbpas
                ti1 = tpg1+(tpg2-tpg1)*dble(i-1)/dble(nbpas)
                ti2 = tpg1+(tpg2-tpg1)*dble(i)/dble(nbpas)
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
                                        ti1, dt10, dt21, &
                                        z2, coef_phase, &
                                        dmoins, metaCurr(nbPhase+SIZE_GRAIN))
                if (metaSteelPara%l_grain_size) then
                    zaust = z2
                    dmoins = metaCurr(nbPhase+SIZE_GRAIN)
                end if
            end do
        else
            dt21_mod = dt21
            taux = metaSteelPara%taux_1+(metaSteelPara%taux_3-metaSteelPara%taux_1)*zeq1
            if ((tpg1 .lt. (metaSteelPara%ac1-epsi)) .or. (zaust .ge. un)) then
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
                                    tpg1, dt10, dt21, &
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
