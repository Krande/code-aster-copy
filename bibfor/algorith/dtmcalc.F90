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

subroutine dtmcalc(sd_dtm_, sd_int_)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! dtmcalc : Main subroutine for the integration of the dynamic equations of motion
!           of a system in generalized coordinates (reduced -modal- basis).
!
#include "jeveux.h"
#include "asterc/etausr.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dtmarch.h"
#include "asterfort/dtmbuff.h"
#include "asterfort/dtmconc.h"
#include "asterfort/dtmdetect.h"
#include "asterfort/dtmget.h"
#include "asterfort/dtmintg.h"
#include "asterfort/dtmupmat.h"
#include "asterfort/intbackup.h"
#include "asterfort/infniv.h"
#include "asterfort/intbuff.h"
#include "asterfort/intget.h"
#include "asterfort/intsav.h"
#include "asterfort/mdidisvisc.h"
#include "asterfort/mdidisisot.h"
#include "asterfort/mdsize.h"
#include "asterfort/nlget.h"
#include "asterfort/resu74.h"
#include "asterfort/sigusr.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpu.h"
#include "asterfort/uttcpr.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_int_
!
!   -0.2- Local variables
    aster_logical     :: upmat, checkcpu
    integer(kind=8)           :: nbrede, nbrevi, nbpas, nbnli, nbmode
    integer(kind=8)           :: n100, i, append, exgyro, lev
    integer(kind=8)           :: adapt, pasarch, iarch, force_arch, reinteg
    integer(kind=8)           :: nltreat, nr, nbsauv, iarch_sd, i_nbar
    integer(kind=8)           :: perc, last_prperc, freqpr, ifm, niv
    integer(kind=8)           :: oldarch, i_nbarf, nbdvis, nbdecr
    real(kind=8)      :: tinit, dt, tps1(4), rint1, rint2
    real(kind=8)      :: time, lastarch, tfin, epsi, newdt
    real(kind=8)      :: dt0
    character(len=8)  :: sd_dtm, sd_int, calcres, nomres, sd_nl
!
    integer(kind=8), pointer :: isto(:) => null()
    integer(kind=8), pointer :: allocs(:) => null()
    integer(kind=8), pointer :: buffdtm(:) => null()
    integer(kind=8), pointer :: buffint(:) => null()
    real(kind=8), pointer :: archlst(:) => null()
    real(kind=8), pointer :: chosav0(:) => null()
    integer(kind=8), pointer :: buffnl(:) => null()

!
!   0 - Initializations
    call infniv(ifm, niv)
!
    sd_dtm = sd_dtm_
    sd_int = sd_int_
    epsi = 10.d0*r8prem()
    force_arch = 0
    iarch = 1
!
!   1 - Retrieval of the necessary information
    call dtmget(sd_dtm, _NB_STEPS, iscal=nbpas)
    call dtmget(sd_dtm, _INST_INI, rscal=tinit)
    call dtmget(sd_dtm, _INST_FIN, rscal=tfin)
!
    call dtmget(sd_dtm, _AR_LINST, vr=archlst)
    call dtmget(sd_dtm, _ARCH_PER, iscal=pasarch)
    call dtmget(sd_dtm, _ADAPT, iscal=adapt)
    call dtmget(sd_dtm, _IARCH_SD, iscal=iarch_sd)
!
    call intget(sd_int, STEP, iocc=1, rscal=dt0)
    time = tinit
    lastarch = tinit

    call uttcpu('CPU.DTMCALC', 'INIT', ' ')
    n100 = nbpas/100+1
!
!   --- From this point on, we loop over the integration steps, extreme care is to be
!       taken in the following as the number of steps can be very large. Avoid thus
!       unnecessary computations inside this loop
!
!   Before looping find out the following :
!   -> Is it necessary to update the matrices on each step
    upmat = .false.
    call dtmget(sd_dtm, _GYRO_FUL, lonvec=exgyro)
    if (exgyro .gt. 0) upmat = .true.

    call dtmget(sd_dtm, _NL_TREAT, iscal=nltreat)
!
!   -> Should we mesure and check the remaining cpu time
    checkcpu = .false.

    call dtmget(sd_dtm, _ARCH_NB, iscal=nbsauv)

    call dtmget(sd_dtm, _APPND_SD, iscal=append)
    call dtmget(sd_dtm, _RESU_SD, kscal=nomres)

!   Already archived steps :  continue case => total number of existing fields
!                          no continue case => a single save (ord = 0) corresponding to
!                                              the initial step
    oldarch = append

    call dtmbuff(sd_dtm, buffdtm)
    call intbuff(sd_int, buffint, level=1)

!   --- Number of forced irregular archiving
    i_nbarf = 0

    perc = 0
    i_nbar = oldarch
    lastarch = tinit
    call utmess('I', 'PROGRESS_1', ni=2, vali=[perc, i_nbar], &
                nr=2, valr=[time, lastarch])
    last_prperc = perc
!
    i = 1
!   --- 5 percent printing period
    freqpr = 5
    if (niv .eq. 2) freqpr = 1
!
    dt = dt0
!   --- Loop until the simulated time reaches tfin with a cumulated numerical error
!       proportional to the number of steps
    do while ((tfin-time) .gt. (i*epsi))
!
        perc = int(100.d0*(time-tinit)/(tfin-tinit))
        if (perc .ne. last_prperc) then
            if (mod(perc, freqpr) .eq. 0) then
                call utmess('I', 'PROGRESS_1', ni=2, vali=[perc, i_nbar], &
                            nr=2, valr=[time, lastarch])
                last_prperc = perc
            end if
        end if
!
!
        if (((i_nbar-oldarch+1) .eq. nbsauv) .and. (iarch_sd .eq. 0)) goto 30

        if (checkcpu) then
!           --- Start recording the integration time required per step
            if (mod(i, n100) .eq. 0) call uttcpu('CPU.DTMCALC', 'DEBUT', ' ')
        end if
!
!       -- Matrices update (Gyroscopy)
        if (upmat) call dtmupmat(sd_dtm, sd_int, buffdtm, buffint)

        if (nltreat .eq. 1) then
            call intbackup(sd_int, '&&INTBAK')
            call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl, buffer=buffdtm)
            call dtmget(sd_dtm, _NL_BUFFER, vi=buffnl, buffer=buffdtm)
            call dtmget(sd_dtm, _NL_SAVE0, vr=chosav0, buffer=buffdtm)
            call nlget(sd_nl, _INTERNAL_VARS, rvect=chosav0, buffer=buffnl)
        end if
!
        nr = 0
10      continue
        nr = nr+1
        if (nr .gt. 100) then
            ASSERT(.false.)
        end if

        call intget(sd_int, STEP, iocc=1, rscal=dt)
        if ((time+dt-archlst(iarch)) .gt. (i*epsi)) then
            newdt = archlst(iarch)-time
            call intsav(sd_int, STEP, 1, iocc=1, rscal=newdt)
            force_arch = 1
        else if (abs(time+dt-archlst(iarch)) .lt. (i*epsi)) then
            newdt = dt
            force_arch = 1
        end if

!       --- Actual integration scheme
        call dtmintg(sd_dtm, sd_int, buffdtm, buffint)
!
!       --- Non-linearity special treatment, is a new reintegration needed
!           (case of a change in state)
        reinteg = 0
        if (nltreat .eq. 1) call dtmdetect(sd_dtm, sd_int, buffdtm, buffint, reinteg)
        if (reinteg .eq. 1) goto 10

        if (i .eq. 1) call intget(sd_int, IND_ARCH, iscal=lev)
        call intget(sd_int, TIME, iocc=lev, rscal=time, buffer=buffint)

        if (force_arch .eq. 2) then
            force_arch = 0
            if (adapt .eq. 0) then
                call intsav(sd_int, STEP, 1, iocc=1, rscal=dt0)
                force_arch = 0
            end if
        end if

        if (force_arch .eq. 1) then
            if (abs(time-archlst(iarch)) .le. (i*epsi)) then
                call dtmarch(sd_dtm, sd_int, buffdtm, buffint)
                iarch = iarch+1
                i_nbar = i_nbar+1
                lastarch = time
                force_arch = 0
                if (abs(dt-newdt) .gt. (i*epsi)) then
                    i_nbarf = i_nbarf+1
                    call intsav(sd_int, STEP, 1, iocc=1, rscal=dt-newdt)
                    force_arch = 2
                else if (adapt .eq. 0) then
                    call intsav(sd_int, STEP, 1, iocc=1, rscal=dt0)
                end if
            end if
        else
            if (mod(i-i_nbarf, pasarch) .eq. 0) then
                i_nbar = i_nbar+1
                call dtmarch(sd_dtm, sd_int, buffdtm, buffint)
                lastarch = time
            end if
        end if

        i = i+1

        if (nr .gt. 1) call intsav(sd_int, STEP, 1, iocc=1, rscal=dt0)

!       --- Interruption by signal usr1 ---
        if (etausr() .eq. 1) call sigusr()

        if (checkcpu) then
!           --- Stop recording the per-step integration time
            if (mod(i, n100) .eq. 0) then
                call uttcpu('CPU.DTMCALC', 'FIN', ' ')
                call uttcpr('CPU.DTMCALC', 4, tps1)
!               --- 1) tps1(4) is the elapsed time (for 1 step)
!                   2) tps1(1) is the remaining time
                rint1 = 5.d0
                rint2 = 0.90d0
!               --- The remaining cpu-time should allow to integrate a chunk of 1% of the
!                   integration steps, if not, resize and finalize the result data structure
!                   then raise an exception
                if (max(rint1, n100*tps1(4)) .gt. (rint2*tps1(1))) then
                    call dtmget(sd_dtm, _CALC_SD, kscal=calcres, buffer=buffdtm)
                    call dtmget(sd_dtm, _NB_MODES, iscal=nbmode, buffer=buffdtm)
                    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnli, buffer=buffdtm)
                    if (nbnli .gt. 0) then
                        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl)
                        call nlget(sd_nl, _NB_REL_FX, iscal=nbrede)
                        call nlget(sd_nl, _NB_REL_FX, iscal=nbrevi)
                    end if
                    call dtmget(sd_dtm, _ARCH_STO, vi=isto, buffer=buffdtm)
!                   --- resize according to the last archived step, isto(1)
                    call mdsize(calcres, isto(1), nbmode, nbnli)
!                   --- Concatenate results in the case of an adaptative integration scheme
                    ! if (adapt.gt.0) then
                    !     call dtmconc(sd_dtm)
                    ! endif
!                   --- Append to an existing result in the case of "reuse"
                    if (append .gt. 0) then
                        call resu74(nomres, calcres)
                    end if

!                   --- Alert the user, raise exceptions and stop the calculation
                    call utmess('Z', 'ALGORITH16_77', ni=2, vali=[i, isto(1)], nr=3, &
                                valr=[lastarch, tps1(4), tps1(1)], num_except=ASTER_TIMELIMIT_ERROR)
                    call utmess('F', 'ALGORITH5_24')
                end if
            end if
        end if
    end do

30  continue
!
    if (last_prperc .ne. 100) then
        perc = 100
        call utmess('I', 'PROGRESS_1', ni=2, vali=[perc, i_nbar], &
                    nr=2, valr=[tfin, lastarch])
    end if

!   --- Concatenate the results if needed
    if (iarch_sd .gt. 0) call dtmconc(sd_dtm)

!   --- Append to an existing result in the case of "reuse"
    if (append .gt. 0) then
        call dtmget(sd_dtm, _CALC_SD, kscal=calcres)
        call resu74(nomres, calcres)
    end if

!   --- Cleanup extra objects
    if (nltreat .eq. 1) then
        call detrsd(' ', '&&DTMMOD')
        call detrsd(' ', '&&DTMNUG')
    end if

    call dtmget(sd_dtm, _NB_NONLI, iscal=nbnli)
    if (nbnli .gt. 0) then
        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl)

        call nlget(sd_nl, _NB_DIS_VISC, iscal=nbdvis)
        call nlget(sd_nl, _NB_DIS_ECRO_TRAC, iscal=nbdecr)
        if ((nbdvis+nbdecr) .gt. 0) then
            call dtmget(sd_dtm, _IND_ALOC, vi=allocs)
            if (nbdvis .gt. 0) then
                call mdidisvisc(sd_nl, nbnli, nomres)
            end if
            if (nbdecr .gt. 0) then
                call mdidisisot(sd_nl, nbnli, nomres)
            end if
        end if

        call detrsd(' ', sd_nl)
    end if
!

end subroutine
