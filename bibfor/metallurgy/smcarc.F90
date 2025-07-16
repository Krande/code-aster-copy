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
subroutine smcarc(nb_hist, nbPhase, ftrc, trc, &
                  coef, fmod, &
                  metaSteelPara, &
                  tempCurr, tPoint, deltaTime, &
                  metaPrev, metaCurr)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterfort/smcaba.h"
#include "asterfort/smcavo.h"
#include "asterfort/smcomo.h"
#include "asterfort/metaSteelTRCPolynom.h"
#include "asterfort/metaSteelGrainSize.h"
#include "asterfort/Metallurgy_type.h"
!
    integer(kind=8), intent(in) :: nb_hist, nbPhase
    real(kind=8), intent(inout) :: ftrc((3*nb_hist), 3), trc((3*nb_hist), 5)
    real(kind=8), intent(in)  :: coef(*), fmod(*)
    type(META_SteelParameters), intent(in) :: metaSteelPara
    real(kind=8), intent(in) :: tempCurr, tPoint, deltaTime
    real(kind=8), intent(in) :: metaPrev(:)
    real(kind=8), intent(out) :: metaCurr(:)
!
! --------------------------------------------------------------------------------------------------
!
! METALLURGY -  Compute phase (steel)
!
! Compute phases (colding)
!
! --------------------------------------------------------------------------------------------------
!
! In  nb_hist             : number of graph in TRC diagram
! In  nbPhase             : number of phases
! IO  trc                 : values of functions for TRC diagram
! IO  ftrc                : values of derivatives (by temperature) of functions for TRC diagram
! In  coef                : parameters from TRC diagrams (P5 polynom)
! In  fmod                : experimental function from TRC diagrams
! In  metaSteelPara       : parameters for metallurgy of steel
! In  tempCurr            : current temperature
! In  tPoint              : increment of temperature
! In  deltaTime           : increment of time
! In  metaPrev            : internal state variables at begin of time step
! Out metaCurr            : internal state variables at end of time step
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0, epsi = 1.d-10
    integer(kind=8) :: ind(6)
    real(kind=8) :: tmf, x(5), tpli
    real(kind=8) :: temp_incr_eff
    real(kind=8) :: coef_phase, incrTempMart
    real(kind=8) :: sumColdPrev, zaust, sumColdCurr, sumFerritePrev, dz(4), dzcold, zMartPrev
    real(kind=8) :: tplm, austeniteMin, akm, bkm, dref, a, zAustPrev
!
! --------------------------------------------------------------------------------------------------
!
    austeniteMin = metaSteelPara%trc%martensiteLaw%austeniteMin
    akm = metaSteelPara%trc%martensiteLaw%akm
    bkm = metaSteelPara%trc%martensiteLaw%bkm
    tplm = metaSteelPara%trc%martensiteLaw%lowerSpeed
    dref = metaSteelPara%trc%austeniteGrain%dref
    a = metaSteelPara%trc%austeniteGrain%a

! - Previous sum of cold phases
    sumColdPrev = metaPrev(PFERRITE)+metaPrev(PPERLITE)+ &
                  metaPrev(PBAINITE)+metaPrev(PMARTENS)
    zAustPrev = 1.d0-sumColdPrev
!
    if (tempCurr .gt. metaSteelPara%ar3) then
! ----- Nothing changes
        metaCurr(PFERRITE) = metaPrev(PFERRITE)
        metaCurr(PPERLITE) = metaPrev(PPERLITE)
        metaCurr(PBAINITE) = metaPrev(PBAINITE)
        metaCurr(PMARTENS) = metaPrev(PMARTENS)
        metaCurr(nbPhase+TEMP_MARTENSITE) = metaSteelPara%ms0
    else

! ----- Température de transformation martensitique courante Ms(t)
        tmf = metaPrev(nbPhase+TEMP_MARTENSITE)-(log(0.01d0))/metaSteelPara%alpha-15.d0

        if ((sumColdPrev .ge. 0.999d0) .or. (tempCurr .lt. tmf)) then
! --------- Nothing changes: hot phase only
            metaCurr(PFERRITE) = metaPrev(PFERRITE)
            metaCurr(PPERLITE) = metaPrev(PPERLITE)
            metaCurr(PBAINITE) = metaPrev(PBAINITE)
            metaCurr(PMARTENS) = metaPrev(PMARTENS)
            metaCurr(nbPhase+TEMP_MARTENSITE) = metaPrev(nbPhase+TEMP_MARTENSITE)
        else
! --------- Compute increment of phases
            if (tempCurr .lt. metaPrev(nbPhase+TEMP_MARTENSITE)) then
                dz(PFERRITE) = zero
                dz(PPERLITE) = zero
                dz(PBAINITE) = zero
            else
! ------------- Compute Teff (effective cooling speed of temperature)
                if (a .eq. 0.d0) then
                    temp_incr_eff = tPoint
                else
                    temp_incr_eff = tPoint*exp(a*(metaPrev(nbPhase+SIZE_GRAIN)-dref))
                end if

! ------------- Compute functions from TRC diagram
                call smcomo(coef, fmod, tempCurr, nb_hist, &
                            ftrc, trc)

                if (temp_incr_eff .gt. (trc(1, 4)*(un+epsi))) then
! ----------------- Before first value from TRC diagrams
                    dz(PFERRITE) = ftrc(1, PFERRITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)
                    dz(PPERLITE) = ftrc(1, PPERLITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)
                    dz(PBAINITE) = ftrc(1, PBAINITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)

                elseif (temp_incr_eff .lt. (trc(nb_hist, 4)*(un-epsi))) then
! ----------------- After last value from TRC diagrams
                    dz(PFERRITE) = &
                        ftrc(nb_hist, PFERRITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)
                    dz(PPERLITE) = &
                        ftrc(nb_hist, PPERLITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)
                    dz(PBAINITE) = &
                        ftrc(nb_hist, PBAINITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)
                else
! ----------------- Find the six nearest TRC curves
                    x(1) = metaPrev(PFERRITE)
                    x(2) = metaPrev(PPERLITE)
                    x(3) = metaPrev(PBAINITE)
                    x(4) = temp_incr_eff
                    x(5) = tempCurr
                    call smcavo(x, nb_hist, trc, ind)

! ----------------- Compute barycenter and update increments of phases
                    call smcaba(x, nb_hist, trc, ftrc, ind, &
                                dz)
                    if ((metaCurr(nbPhase+STEEL_TEMP)-tempCurr) .gt. zero) then
                        dz(PFERRITE) = zero
                        dz(PPERLITE) = zero
                        dz(PBAINITE) = zero
                    else
                        dz(PFERRITE) = dz(PFERRITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)
                        dz(PPERLITE) = dz(PPERLITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)
                        dz(PBAINITE) = dz(PBAINITE)*(metaCurr(nbPhase+STEEL_TEMP)-tempCurr)
                    end if
                end if
            end if

! --------- New value of sum of three first phases
            sumFerritePrev = sumColdPrev-metaPrev(PMARTENS)

! --------- Compute new martensite temperature
            if ((sumFerritePrev .ge. austeniteMin) .and. (metaPrev(PMARTENS) .eq. zero)) then
                metaCurr(nbPhase+TEMP_MARTENSITE) = metaSteelPara%ms0+akm*sumFerritePrev+bkm
            else
                metaCurr(nbPhase+TEMP_MARTENSITE) = metaPrev(nbPhase+TEMP_MARTENSITE)
            end if

! --------- Compute proportion of martensite
            zMartPrev = 1.d0-sumFerritePrev
            if ((metaCurr(nbPhase+STEEL_TEMP) .gt. metaCurr(nbPhase+TEMP_MARTENSITE)) .or. &
                (zMartPrev .lt. 0.01d0)) then
                metaCurr(PMARTENS) = metaPrev(PMARTENS)
            else
! ------------- Compute derivative of temperature
                call metaSteelTRCPolynom(coef(3:8), tplm, tempCurr, tpli)
                if ((tPoint .gt. tpli) .and. (metaPrev(PMARTENS) .eq. zero)) then
                    metaCurr(PMARTENS) = metaPrev(PMARTENS)
                else

                    incrTempMart = &
                        max(0.d0, metaCurr(nbPhase+TEMP_MARTENSITE)-metaCurr(nbPhase+STEEL_TEMP))

                    metaCurr(PMARTENS) = zMartPrev* &
                                         (un-exp(metaSteelPara%alpha*incrTempMart))
                end if
            end if
            dz(PMARTENS) = metaCurr(PMARTENS)-metaPrev(PMARTENS)

! --------- Always create martensite
            if (dz(PMARTENS) .lt. 0.d0) then
                dz(PMARTENS) = 0.d0
            end if

! --------- New value of sum of "cold" phases
            sumColdCurr = sumColdPrev+dz(PFERRITE)+dz(PPERLITE)+ &
                          dz(PBAINITE)+dz(PMARTENS)

            if (sumColdCurr .gt. 0.999d0) then
                dzcold = sumColdCurr-sumColdPrev
                metaCurr(PFERRITE) = metaPrev(PFERRITE)+dz(PFERRITE)/(dzcold/(un-sumColdPrev))
                metaCurr(PPERLITE) = metaPrev(PPERLITE)+dz(PPERLITE)/(dzcold/(un-sumColdPrev))
                metaCurr(PBAINITE) = metaPrev(PBAINITE)+dz(PBAINITE)/(dzcold/(un-sumColdPrev))
                metaCurr(PMARTENS) = metaPrev(PMARTENS)+dz(PMARTENS)/(dzcold/(un-sumColdPrev))
            else
                metaCurr(PFERRITE) = metaPrev(PFERRITE)+dz(PFERRITE)
                metaCurr(PPERLITE) = metaPrev(PPERLITE)+dz(PPERLITE)
                metaCurr(PBAINITE) = metaPrev(PBAINITE)+dz(PBAINITE)
                metaCurr(PMARTENS) = metaPrev(PMARTENS)+dz(PMARTENS)
            end if
        end if
    end if

! - Compute "hot" phase
    zaust = un-metaCurr(PFERRITE)+metaCurr(PPERLITE)+ &
            metaCurr(PBAINITE)+metaCurr(PMARTENS)

! - Compute size of grain
    coef_phase = un
    call metaSteelGrainSize(metaSteelPara, &
                            tempCurr, deltaTime, deltaTime, &
                            zaust, coef_phase, &
                            metaPrev(nbPhase+SIZE_GRAIN), metaCurr(nbPhase+SIZE_GRAIN))
!
end subroutine
