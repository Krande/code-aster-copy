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
subroutine intadapt2(sd_dtm_, sd_int_, buffdtm, buffint)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! intadapt2 : Integrate from t_i to t_i+1 the differential equations of motion
!             using the ADAPT order-2 formulation.
!
#include "jeveux.h"
#include "blas/dcopy.h"
#include "asterc/r8prem.h"
#include "asterfort/dtmacce.h"
#include "asterfort/dtmget.h"
#include "asterfort/frqapp.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/intbuff.h"
#include "asterfort/intget.h"
#include "asterfort/intinivec.h"
#include "asterfort/intsav.h"
#include "asterfort/nlget.h"
#include "asterfort/utmess.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_int_
    integer(kind=8), pointer :: buffdtm(:)
    integer(kind=8), pointer :: buffint(:)
!
!   -0.2- Local variables
    integer(kind=8) :: i, nbequ, ind1, iret1, iret2
    integer(kind=8) :: npper, nrmax, nr, nbnoli, nbvint
    real(kind=8) :: t1, dt, dt1, dt2, pas1
    real(kind=8) :: pas2, coeff, err, epsi, norm
    real(kind=8) :: cpmin, freq
    character(len=8) :: sd_dtm, sd_int, sd_nl, vvar
!
    real(kind=8), pointer :: depl1(:) => null()
    real(kind=8), pointer :: vite1(:) => null()
    real(kind=8), pointer :: acce1(:) => null()
    real(kind=8), pointer :: fext1(:) => null()
    real(kind=8), pointer :: depl2(:) => null()
    real(kind=8), pointer :: vite2(:) => null()
    real(kind=8), pointer :: acce2(:) => null()
    real(kind=8), pointer :: fext2(:) => null()
!
    real(kind=8), pointer :: par(:) => null()
    real(kind=8), pointer :: vmin(:) => null()
    real(kind=8), pointer :: velint(:) => null()
!
    real(kind=8), pointer :: nlsav0(:) => null()
    real(kind=8), pointer :: nlsav1(:) => null()
!
    integer(kind=8), pointer :: buffnl(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! Algorithm parameters saved in a linear vector, easily accessible using the defines
#define crit_vmin par(1)
#define npper_r par(2)
#define nrmax_r par(3)
#define cmult par(4)
#define cdivi par(5)
#define dtmin par(6)
#define dtmax par(7)
#define deltadt par(8)
#define stabstep par(9)
#define nbnlsav par(10)
#define nbsavnl nint(par(10))
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_int = sd_int_
    epsi = 100.d0*r8prem()
!
!   1 - Retrieval of the system's state at instant t_i (index=1)
    call intget(sd_int, TIME, iocc=1, rscal=t1, buffer=buffint)
    call intget(sd_int, INDEX, iocc=1, iscal=ind1, buffer=buffint)
    call intget(sd_int, STEP, iocc=1, rscal=dt, buffer=buffint)
    dt1 = dt
!
    call intget(sd_int, DEPL, iocc=1, vr=depl1, lonvec=nbequ, &
                buffer=buffint)
    call intget(sd_int, VITE, iocc=1, vr=vite1, buffer=buffint)
    call intget(sd_int, ACCE, iocc=1, vr=acce1, buffer=buffint)
    call intget(sd_int, FORCE_EX, iocc=1, vr=fext1, buffer=buffint)
!
!   2 - Detection of the initial call to the ADAPT_ORDR2 algorithm
!       DEPL/2 does not exist in the buffer
    call intget(sd_int, DEPL, iocc=2, lonvec=iret1, buffer=buffint)
    if (iret1 .eq. 0) then
!
!       2.1 - Algorithm initialization
!
        call intinivec(sd_int, PARAMS, 10, vr=par)
!
        call getvtx('SCHEMA_TEMPS', 'VITE_MIN', iocc=1, scal=vvar)
!       VITE_MIN = 'NORM' <=> crit_vmin = 1.d0 ; VITE_MIN = 'MAXI' <=> crit_vmin = 2.d0
        crit_vmin = 1.d0
        if (vvar(1:4) .eq. 'MAXI') crit_vmin = 2.d0
!
        call getvis('SCHEMA_TEMPS', 'NB_POIN_PERIODE', iocc=1, scal=npper)
        call getvis('SCHEMA_TEMPS', 'NMAX_ITER_PAS', iocc=1, scal=nrmax)
        npper_r = 1.d0*npper
        nrmax_r = 1.d0*nrmax
!
        call getvr8('SCHEMA_TEMPS', 'COEF_MULT_PAS', iocc=1, scal=cmult)
        call getvr8('SCHEMA_TEMPS', 'COEF_DIVI_PAS', iocc=1, scal=cdivi)
!
        call getvr8('SCHEMA_TEMPS', 'PAS_MAXI', iocc=1, scal=dtmax, nbret=iret1)
        if (iret1 .ne. 1) dtmax = 1.d6*dt
        call getvr8('SCHEMA_TEMPS', 'PAS_MINI', iocc=1, scal=dtmin, nbret=iret2)
        if (iret2 .ne. 1) then
            call getvr8('SCHEMA_TEMPS', 'PAS_LIMI_RELA', iocc=1, scal=cpmin)
            dtmin = cpmin*dt
        end if
!
!       deltadt gives the ratio between dtmin and dtmax, it is considered as an
!       indicator for whether we should adapt or no the time step
        deltadt = abs(dtmax/dtmin-1.d0)
!
        stabstep = 1.d0
!
!       --- Allocate and initialize work vector vmin
        call intinivec(sd_int, WORK1, nbequ, vr=vmin)
        do i = 1, nbequ
            vmin(i) = 1.d-15
        end do
!
!       --- Allocate work vector giving intermediate velocity at t_i+1/2
        call intinivec(sd_int, WORK2, nbequ, vr=velint)
!
!       --- Allocate work vectors for NL_SAVES
        call dtmget(sd_dtm, _NB_NONLI, iscal=nbnoli, buffer=buffdtm)
        if (nbnoli .gt. 0) then
            call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl, buffer=buffdtm)
            call nlget(sd_nl, _INTERNAL_VARS, lonvec=nbvint)
            nbnlsav = nbvint*1.d0
            call intinivec(sd_int, WORK3, nbvint, vr=nlsav0)
        else
            nbnlsav = 0.d0
        end if
!
!       --- Allocate vectors DEPL/VITE/ACCE/2 (t_i+1)
        call intinivec(sd_int, DEPL, nbequ, iocc=2, vr=depl2)
        call intinivec(sd_int, VITE, nbequ, iocc=2, vr=vite2)
        call intinivec(sd_int, ACCE, nbequ, iocc=2, vr=acce2)
        call intinivec(sd_int, FORCE_EX, nbequ, iocc=2, vr=fext2)
!
        dt = 0.d0
!
        nullify (buffint)
        call intbuff(sd_int, buffint, level=2)
    else
!       --- Algorithm is already initialized, just retrieve DEPL/VITE/ACCE/2
!           Note that these are made identical to index 1 (dcopy before exit)
        call intget(sd_int, DEPL, iocc=2, vr=depl2, buffer=buffint)
        call intget(sd_int, VITE, iocc=2, vr=vite2, buffer=buffint)
        call intget(sd_int, ACCE, iocc=2, vr=acce2, buffer=buffint)
!
!       --- Retrieve algorithm parameters
        call intget(sd_int, PARAMS, vr=par, buffer=buffint)
!       --- Retrieve work vectors vmin, and vinter
        call intget(sd_int, WORK1, vr=vmin, buffer=buffint)
        call intget(sd_int, WORK2, vr=velint, buffer=buffint)
!       --- Retrieve choc parameters save container
        if (nbsavnl .gt. 0) call intget(sd_int, WORK3, vr=nlsav0, buffer=buffint)
    end if
    if (nbsavnl .gt. 0) then
        call dtmget(sd_dtm, _SD_NONL, kscal=sd_nl, buffer=buffdtm)
        call dtmget(sd_dtm, _NL_BUFFER, vi=buffnl, buffer=buffdtm)
        call nlget(sd_nl, _INTERNAL_VARS, vr=nlsav1, buffer=buffnl)
        b_n = to_blas_int(nbsavnl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, nlsav1, b_incx, nlsav0, b_incy)
    end if
    call intsav(sd_int, INDEX, 1, iocc=2, iscal=ind1+1, &
                buffer=buffint)
!
    nr = 1
!
10  continue
!
    if (nbsavnl .gt. 0) then
        b_n = to_blas_int(nbsavnl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, nlsav0, b_incx, nlsav1, b_incy)
    end if
!
!   3 - Calculate the system's state at index <2>
!
    pas1 = 0.5d0*(dt+dt1)
    pas2 = 0.5d0*dt1
!
    do i = 1, nbequ
        velint(i) = vite1(i)+acce1(i)*pas1
        depl2(i) = depl1(i)+velint(i)*dt1
        vite2(i) = velint(i)+pas2*acce1(i)
    end do
!
    call intsav(sd_int, TIME, 1, iocc=2, rscal=t1+dt1, &
                buffer=buffint)
    call intsav(sd_int, STEP, 1, iocc=2, rscal=dt1, &
                buffer=buffint)
!
    call dtmacce(sd_dtm, sd_int, 2, buffdtm, buffint)
!
!   4 - If an adaptative scheme is requested (dtmax != dtmin)
    if (deltadt .ge. epsi) then
!
!       4.1 - Estimation of the error
        call frqapp(dt1, nbequ, depl1, depl2, acce1, &
                    acce2, vmin, freq)
!
        err = max(npper_r*freq*dt1, epsi)
!
        coeff = 1.d0
        if (err .ge. 1.d0) then
            coeff = 1.d0/cdivi
        else if (stabstep .gt. 4.5d0) then
            if (err .lt. 0.75d0) then
                coeff = cmult
                stabstep = 1.d0
            end if
        else
            stabstep = stabstep+1.d0
        end if
!
!       4.2 - Determine the time step for the next iteration or integration step
        dt2 = min(dtmax, dt1*coeff)
!
        if ((err .ge. 1.d0) .and. (dt2 .ge. dtmin) .and. (nr .lt. nint(nrmax_r))) then
            nr = nr+1
            dt1 = dt2
            goto 10
        else if (dt2 .lt. dtmin) then
            call utmess('F', 'ALGORITH5_23')
        else if (nr .eq. nint(nrmax_r)) then
            call utmess('A', 'DYNAMIQUE_18', si=nr, nr=2, valr=[t1, dt1])
        end if
!
!       4.3 - Caculate Vmin for the next step
        if (abs(crit_vmin-1.d0) .le. epsi) then
!           --- NORM
            norm = 0.d0
            do i = 1, nbequ
                norm = norm+velint(i)**2
            end do
            norm = sqrt(norm)*0.01d0
            do i = 1, nbequ
                vmin(i) = norm
            end do
        else
!           --- MAXI
            do i = 1, nbequ
                vmin(i) = max(abs(velint(i)*0.01d0), vmin(i))
            end do
        end if
!
    end if
!
!   5 - Preparing the algorithm for the next step, copy index 2 in 1
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, depl2, b_incx, depl1, b_incy)
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, velint, b_incx, vite1, b_incy)
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, acce2, b_incx, acce1, b_incy)
    call intsav(sd_int, STEP, 1, iocc=1, rscal=dt2, &
                buffer=buffint)
    call intsav(sd_int, TIME, 1, iocc=1, rscal=t1+dt1, &
                buffer=buffint)
!
    call intsav(sd_int, INDEX, 1, iocc=1, iscal=ind1+1, &
                buffer=buffint)
    call intsav(sd_int, INDEX, 1, iocc=2, iscal=ind1+1, &
                buffer=buffint)
!
!   6 - Set the archiving index to 2
    call intsav(sd_int, IND_ARCH, 1, iscal=2, buffer=buffint)
!
end subroutine
