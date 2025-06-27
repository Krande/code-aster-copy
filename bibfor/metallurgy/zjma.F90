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
subroutine zjma(metaSteelPara, &
                nbVari, nbVariTemper, nbVariPrev, &
                temp1, temp2, &
                deltaTime12, &
                prevMetaIsTemper, &
                metaPrev, metaCurr, metaCurrTemper)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
!
    integer(kind=8), intent(in) :: nbVari, nbVariTemper, nbVariPrev
    type(META_SteelParameters), intent(in) :: metaSteelPara
    real(kind=8), intent(in) :: deltaTime12, temp1, temp2
    aster_logical, intent(in) :: prevMetaIsTemper
    real(kind=8), intent(in) :: metaPrev(nbVariPrev)
    real(kind=8), intent(in) :: metaCurr(nbVari)
    real(kind=8), intent(out) :: metaCurrTemper(nbVariTemper)
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY -  Tempering law for steel (bainite and martensite) phase computing
!
! --------------------------------------------------------------------------------------------------
!
! In  metaSteelPara      : parameters for metallurgy of steel
! In  metaPrev           : value of internal state variable at previous time step
! Out metaCurr           : value of internal state variable at current time step
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: tempPgPrev
    real(kind=8) :: deltaTemp
    integer(kind=8) :: cyclTherPrev, cyclTherCurr
    real(kind=8) :: ZTildeMartPrev, ZTildeBainPrev
    real(kind=8) :: ZTildeMartCurr, ZTildeBainCurr
    real(kind=8) :: tau_0_bain, tau_0_mart
    real(kind=8) :: tempTempering, tempPgCurr
    real(kind=8) :: ZBainBrut, ZMartBrut
    real(kind=8) :: ZBainRevePrev, ZMartRevePrev
    real(kind=8) :: ZBainBaseCurr, ZMartBaseCurr
    real(kind=8) :: ZBainReveCurr, ZMartReveCurr
    real(kind=8) :: ZBainTotaPrev, ZMartTotaPrev
    real(kind=8) :: bainCoefB, bainCoefN
    real(kind=8) :: martCoefB, martCoefN, tempHold
    real(kind=8) :: R, delta_H, C0, deltaTimeEqui, Pa
    real(kind=8), parameter :: kelvin = 273.0d0
!
! --------------------------------------------------------------------------------------------------
!
    deltaTemp = abs(temp2-temp1)

! - Get previous internal state variables
    if (prevMetaIsTemper) then
        tempPgPrev = metaPrev(STEEL_TEMP+NBPHASESTEELR)
        cyclTherPrev = nint(metaPrev(THER_CYCL+NBPHASESTEELR))
        ZBainRevePrev = metaPrev(PRBAINITER)
        ZMartRevePrev = metaPrev(PRMARTENSR)
    else
        tempPgPrev = metaPrev(STEEL_TEMP+NBPHASESTEEL)
        ZBainRevePrev = 0.d0
        ZMartRevePrev = 0.d0
        cyclTherPrev = 0
    end if

! - Get internal state variables without tempering
    tempPgCurr = metaCurr(STEEL_TEMP+NBPHASESTEEL)
    ZBainBrut = metaCurr(PBAINITE)
    ZMartBrut = metaCurr(PMARTENS)

    deltaTemp = abs(tempPgCurr-tempPgPrev)

! - Compute ratio
    ZBainTotaPrev = ZBainRevePrev+ZBainBrut
    ZMartTotaPrev = ZMartRevePrev+ZMartBrut
    ZTildeMartPrev = 0.d0
    if (abs(ZMartTotaPrev) .ge. r8prem()) then
        ZTildeMartPrev = ZMartRevePrev/ZMartTotaPrev
    end if
    ZTildeBainPrev = 0.d0
    if (abs(ZBainTotaPrev) .ge. r8prem()) then
        ZTildeBainPrev = ZBainRevePrev/ZBainTotaPrev
    end if

! - Get material parameters
    tempTempering = metaSteelPara%temper%temp
    tempHold = metaSteelPara%temper%tempHold
    bainCoefB = metaSteelPara%temper%bainite_b
    bainCoefN = metaSteelPara%temper%bainite_n
    martCoefB = metaSteelPara%temper%martensite_b
    martCoefN = metaSteelPara%temper%martensite_n
    R = 2.0
    delta_H = 100000.0
    C0 = (log(10.d0)*R)/delta_H

! - Computes
    cyclTherCurr = cyclTherPrev
    if (abs(deltaTime12) .le. r8prem()) then
        ZBainBaseCurr = ZBainBrut
        ZMartBaseCurr = ZMartBrut
        ZBainReveCurr = ZBainRevePrev
        ZMartReveCurr = ZMartRevePrev
        cyclTherCurr = cyclTherPrev
    else
        ZTildeMartCurr = 0.d0
        ZTildeBainCurr = 0.d0
        if (tempPgCurr .gt. metaSteelPara%ac3) then
            cyclTherCurr = 0
            ZTildeMartCurr = 0.0
            ZTildeBainCurr = 0.0
        else if (tempPgCurr .gt. tempTempering) then
            if (cyclTherPrev .eq. 2) then
! ------------- Bainite
                tau_0_bain = (-1.0*log(1-ZBainRevePrev)/bainCoefB)**(1.0/bainCoefN)
                Pa = 1.0/(1.0/(tempPgCurr+kelvin)-C0*log(deltaTime12))
                deltaTimeEqui = exp((1.0/C0)*((1.0/(tempHold+kelvin))-(1.0/Pa)))
                ZTildeBainCurr = &
                    1.0-exp(-1.0*bainCoefB*(tau_0_bain+deltaTimeEqui)**bainCoefN)
! ------------- Martensite
                tau_0_mart = (-1.0*log(1-ZMartRevePrev)/martCoefB)**(1.0/martCoefN)
                Pa = 1.0/(1.0/(tempPgCurr+kelvin)-C0*log(deltaTime12))
                deltaTimeEqui = exp((1.0/C0)*((1.0/(tempHold+kelvin))-(1.0/Pa)))
                ZTildeMartCurr = &
                    1.0-exp(-1.0*martCoefB*(tau_0_mart+deltaTimeEqui)**martCoefN)
            end if
        else
            if (cyclTherPrev .eq. 2) then
                ZTildeMartCurr = 0.d0
                ZTildeBainCurr = 0.d0
                if (ZMartBrut .ge. r8prem()) then
                    ZTildeMartCurr = ZMartRevePrev/ZMartBrut
                end if
                if (ZBainBrut .ge. r8prem()) then
                    ZTildeBainCurr = ZBainRevePrev/ZBainBrut
                end if
                cyclTherCurr = cyclTherPrev
            elseif (cyclTherPrev .eq. 1) then
                ZTildeMartCurr = 0.d0
                ZTildeBainCurr = 0.d0
                if (deltaTemp .ge. r8prem()) then
                    cyclTherCurr = 2
                end if
            else if (cyclTherPrev .eq. 0) then
                ZTildeMartCurr = 0.d0
                ZTildeBainCurr = 0.d0
                if (ZBainBrut .ge. r8prem() .or. ZMartBrut .ge. r8prem()) then
                    cyclTherCurr = 1
                end if
            else
                WRITE (6, *) "cyclTherPrev: ", cyclTherPrev
                ASSERT(ASTER_FALSE)
            end if
        end if
        ZMartReveCurr = ZMartBrut*ZTildeMartCurr
        ZBainReveCurr = ZBainBrut*ZTildeBainCurr
        ZMartBaseCurr = ZMartBrut-ZMartReveCurr
        ZBainBaseCurr = ZBainBrut-ZBainReveCurr
    end if

! - Update internal state variables
    metaCurrTemper(PRFERRITE) = metaCurr(PFERRITE)
    metaCurrTemper(PRPERLITE) = metaCurr(PPERLITE)
    metaCurrTemper(PRBAINITEB) = ZBainBaseCurr
    metaCurrTemper(PRMARTENSB) = ZMartBaseCurr
    metaCurrTemper(PRBAINITER) = ZBainReveCurr
    metaCurrTemper(PRMARTENSR) = ZMartReveCurr
    metaCurrTemper(PRAUSTENITE) = metaCurr(PAUSTENITE)
    metaCurrTemper(PRSUMCOLD) = 1.d0-metaCurr(PAUSTENITE)
    metaCurrTemper(SIZE_GRAIN+NBPHASESTEELR) = metaCurr(SIZE_GRAIN+NBPHASESTEEL)
    metaCurrTemper(STEEL_TEMP+NBPHASESTEELR) = metaCurr(STEEL_TEMP+NBPHASESTEEL)
    metaCurrTemper(TEMP_MARTENSITE+NBPHASESTEELR) = metaCurr(TEMP_MARTENSITE+NBPHASESTEEL)
    metaCurrTemper(THER_CYCL+NBPHASESTEELR) = cyclTherCurr
!
end subroutine
