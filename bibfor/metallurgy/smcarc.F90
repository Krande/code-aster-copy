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
subroutine smcarc(nb_hist, nb_phase, ftrc, trc, &
                  coef, fmod, &
                  metaSteelPara, &
                  temp_curr, temp_incr, time_incr, &
                  vari_prev, vari_curr)
!
    use Metallurgy_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/smcaba.h"
#include "asterfort/smcavo.h"
#include "asterfort/smcomo.h"
#include "asterfort/metaSteelTRCPolynom.h"
#include "asterfort/metaSteelGrainSize.h"
#include "asterfort/Metallurgy_type.h"
!
    integer, intent(in) :: nb_hist, nb_phase
    real(kind=8), intent(inout) :: ftrc((3*nb_hist), 3), trc((3*nb_hist), 5)
    real(kind=8), intent(in)  :: coef(*), fmod(*)
    type(META_SteelParameters), intent(in) :: metaSteelPara
    real(kind=8), intent(in) :: temp_curr, temp_incr, time_incr
    real(kind=8), intent(in) :: vari_prev(:)
    real(kind=8), intent(out) :: vari_curr(:)
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
! In  nb_phase            : number of phases
! IO  trc                 : values of functions for TRC diagram
! IO  ftrc                : values of derivatives (by temperature) of functions for TRC diagram
! In  coef                : parameters from TRC diagrams (P5 polynom)
! In  fmod                : experimental function from TRC diagrams
! In  metaSteelPara       : parameters for metallurgy of steel
! In  temp_curr           : current temperature
! In  temp_incr           : increment of temperature
! In  time_incr           : increment of time
! In  vari_prev           : internal state variables at begin of time step
! Out vari_curr           : internal state variables at end of time step
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0, epsi = 1.d-10
    integer :: ind(6)
    real(kind=8) :: tmf, x(5), tpli
    real(kind=8) :: temp_incr_eff
    real(kind=8) :: coef_phase
    real(kind=8) :: zcold0, zaust, zcold, z13, dz(4), dzcold, zmartensite
    real(kind=8) :: tplm, austeniteMin, akm, bkm, dref, a
!
! --------------------------------------------------------------------------------------------------
!
    austeniteMin = metaSteelPara%trc%martensiteLaw%austeniteMin
    akm = metaSteelPara%trc%martensiteLaw%akm
    bkm = metaSteelPara%trc%martensiteLaw%bkm
    tplm = metaSteelPara%trc%martensiteLaw%lowerSpeed
    dref = metaSteelPara%trc%austeniteGrain%dref
    a = metaSteelPara%trc%austeniteGrain%a

!
    if (temp_curr .gt. metaSteelPara%ar3) then
! ----- Nothing changes
        vari_curr(PFERRITE) = vari_prev(PFERRITE)
        vari_curr(PPERLITE) = vari_prev(PPERLITE)
        vari_curr(PBAINITE) = vari_prev(PBAINITE)
        vari_curr(PMARTENS) = vari_prev(PMARTENS)
        vari_curr(nb_phase+TEMP_MARTENSITE) = metaSteelPara%ms0
    else
        zcold0 = vari_prev(PFERRITE)+vari_prev(PPERLITE)+ &
                 vari_prev(PBAINITE)+vari_prev(PMARTENS)
        tmf = vari_prev(nb_phase+TEMP_MARTENSITE)-(log(0.01d0))/metaSteelPara%alpha-15.d0
        if ((zcold0 .ge. 0.999d0) .or. (temp_curr .lt. tmf)) then
! --------- Nothing changes: hot phase only
            vari_curr(PFERRITE) = vari_prev(PFERRITE)
            vari_curr(PPERLITE) = vari_prev(PPERLITE)
            vari_curr(PBAINITE) = vari_prev(PBAINITE)
            vari_curr(PMARTENS) = vari_prev(PMARTENS)
            vari_curr(nb_phase+TEMP_MARTENSITE) = vari_prev(nb_phase+TEMP_MARTENSITE)
        else
! --------- Compute increment of phases
            if (temp_curr .lt. vari_prev(nb_phase+TEMP_MARTENSITE)) then
                dz(PFERRITE) = zero
                dz(PPERLITE) = zero
                dz(PBAINITE) = zero
            else
! ------------- Compute Teff (effective cooling speed of temperature)
                if (a .eq. 0.d0) then
                    temp_incr_eff = temp_incr
                else
                    temp_incr_eff = temp_incr*exp(a*(vari_prev(nb_phase+SIZE_GRAIN)-dref))
                end if
! ------------- Compute functions from TRC diagram
                call smcomo(coef, fmod, temp_curr, nb_hist, &
                            ftrc, trc)
                if (temp_incr_eff .gt. (trc(1, 4)*(un+epsi))) then
! ----------------- Before first value from TRC diagrams
                    dz(PFERRITE) = ftrc(1, PFERRITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                    dz(PPERLITE) = ftrc(1, PPERLITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                    dz(PBAINITE) = ftrc(1, PBAINITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                elseif (temp_incr_eff .lt. (trc(nb_hist, 4)*(un-epsi))) then
! ----------------- After last value from TRC diagrams
                    dz(PFERRITE) = &
                        ftrc(nb_hist, PFERRITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                    dz(PPERLITE) = &
                        ftrc(nb_hist, PPERLITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                    dz(PBAINITE) = &
                        ftrc(nb_hist, PBAINITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                else
! ----------------- Find the six nearest TRC curves
                    x(1) = vari_prev(PFERRITE)
                    x(2) = vari_prev(PPERLITE)
                    x(3) = vari_prev(PBAINITE)
                    x(4) = temp_incr_eff
                    x(5) = temp_curr
                    call smcavo(x, nb_hist, trc, ind)
! ----------------- Compute barycenter and update increments of phases
                    call smcaba(x, nb_hist, trc, ftrc, ind, &
                                dz)
                    if ((vari_curr(nb_phase+STEEL_TEMP)-temp_curr) .gt. zero) then
                        dz(PFERRITE) = zero
                        dz(PPERLITE) = zero
                        dz(PBAINITE) = zero
                    else
                        dz(PFERRITE) = dz(PFERRITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                        dz(PPERLITE) = dz(PPERLITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                        dz(PBAINITE) = dz(PBAINITE)*(vari_curr(nb_phase+STEEL_TEMP)-temp_curr)
                    end if
                end if
            end if
! --------- New value of sum of three first phases
            z13 = zcold0-vari_prev(PMARTENS)
! --------- Compute new martensite temperature
            if ((z13 .ge. austeniteMin) .and. (vari_prev(PMARTENS) .eq. zero)) then
                vari_curr(nb_phase+TEMP_MARTENSITE) = metaSteelPara%ms0+akm*z13+bkm
            else
                vari_curr(nb_phase+TEMP_MARTENSITE) = vari_prev(nb_phase+TEMP_MARTENSITE)
            end if
! --------- Compute proportion of martensite
            zmartensite = 1.d0-z13
            if ((vari_curr(nb_phase+STEEL_TEMP) .gt. vari_curr(nb_phase+TEMP_MARTENSITE)) .or. &
                (zmartensite .lt. 0.01d0)) then
                vari_curr(PMARTENS) = vari_prev(PMARTENS)
            else
! ------------- Compute derivative of temperature
                call metaSteelTRCPolynom(coef(3:8), tplm, temp_curr, tpli)
                if ((temp_incr .gt. tpli) .and. (vari_prev(PMARTENS) .eq. zero)) then
                    vari_curr(PMARTENS) = vari_prev(PMARTENS)
                else
                    vari_curr(PMARTENS) = zmartensite* &
                                          (un-exp(metaSteelPara%alpha* &
                                                  (vari_curr(nb_phase+TEMP_MARTENSITE)- &
                                                   vari_curr(nb_phase+STEEL_TEMP))))
                end if
            end if
            dz(PMARTENS) = vari_curr(PMARTENS)-vari_prev(PMARTENS)
! --------- New value of sum of "cold" phases
            zcold = vari_prev(PFERRITE)+vari_prev(PPERLITE)+ &
                    vari_prev(PBAINITE)+vari_prev(PMARTENS)
            zcold = zcold+dz(PFERRITE)+dz(PPERLITE)+ &
                    dz(PBAINITE)+dz(PMARTENS)
            if (zcold .gt. 0.999d0) then
                dzcold = zcold-zcold0
                vari_curr(PFERRITE) = vari_prev(PFERRITE)+dz(PFERRITE)/(dzcold/(un-zcold0))
                vari_curr(PPERLITE) = vari_prev(PPERLITE)+dz(PPERLITE)/(dzcold/(un-zcold0))
                vari_curr(PBAINITE) = vari_prev(PBAINITE)+dz(PBAINITE)/(dzcold/(un-zcold0))
                vari_curr(PMARTENS) = vari_prev(PMARTENS)+dz(PMARTENS)/(dzcold/(un-zcold0))
            else
                vari_curr(PFERRITE) = vari_prev(PFERRITE)+dz(PFERRITE)
                vari_curr(PPERLITE) = vari_prev(PPERLITE)+dz(PPERLITE)
                vari_curr(PBAINITE) = vari_prev(PBAINITE)+dz(PBAINITE)
                vari_curr(PMARTENS) = vari_prev(PMARTENS)+dz(PMARTENS)
            end if
        end if
    end if

    ! if (vari_curr(PFERRITE) .gt. 1.d0) then
    !     WRITE (6, *) "zcold0: ", zcold0
    !     WRITE (6, *) "tmf: ", tmf
    !     WRITE (6, *) "temp_curr: ", temp_curr
    !     WRITE (6, *) "metaSteelPara%ar3: ", metaSteelPara%ar3
    !     WRITE (6, *) "PFERRITe: ", vari_prev(PFERRITE), dz(PFERRITE), vari_curr(PFERRITE)
    !     ASSERT(ASTER_FALSE)
    ! end if

! - Compute "hot" phase
    zaust = un-vari_curr(PFERRITE)+vari_curr(PPERLITE)+ &
            vari_curr(PBAINITE)+vari_prev(PBAINITE)

! - Compute size of grain
    coef_phase = un
    call metaSteelGrainSize(metaSteelPara, &
                            temp_curr, time_incr, time_incr, &
                            zaust, coef_phase, &
                            vari_prev(nb_phase+SIZE_GRAIN), vari_curr(nb_phase+SIZE_GRAIN))
!
end subroutine
