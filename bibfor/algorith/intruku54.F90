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
subroutine intruku54(sd_dtm_, sd_int_, buffdtm, buffint)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! intruku54 : Integrate from t_i to t_i+1 the differential equations of motion
!             using a Runge-Kutta type order 5/4 scheme.
!
#include "jeveux.h"
#include "blas/dcopy.h"
#include "asterc/r8prem.h"
#include "asterc/r8maem.h"
#include "asterfort/dtmacce.h"
#include "asterfort/dtmget.h"
#include "asterfort/getvr8.h"
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
    integer(kind=8) :: ee, ss, j_kde, j_kvi, nbvint
    integer(kind=8) :: nbnoli, iret, iret3
    real(kind=8) :: t1, dt, t2, dt2, coeff
    real(kind=8) :: errd, errde, errv, errvi, errt
    real(kind=8) :: seuil1, seuil2, skd, skv, epsi
    real(kind=8) :: ddep, dvit, errt0
    character(len=8) :: sd_dtm, sd_int, sd_nl
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
    real(kind=8), pointer :: nlsav0(:) => null()
    real(kind=8), pointer :: nlsav1(:) => null()
!
    real(kind=8), pointer :: par(:) => null()
    integer(kind=8), pointer :: buffnl(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
! Algorithm parameters saved in a linear vector, easily accessible using the defines
#define cdp(i) par(i)
#define adp(i,j) par(7+(i-1)*6+j)
#define b4b5(i) par(7+36+i)
#define tol par(7+36+7+1)
#define alpha par(7+36+7+2)
#define dtmin par(7+36+7+3)
#define dtmax par(7+36+7+4)
#define deltadt par(7+36+7+5)
#define nbnlsav par(7+36+7+6)
#define nbsavnl nint(par(7+36+7+6))
! Work vectors accessible using indices (ee,i) where :
! <ee> : RK level (1-6)
! <i>  : equation index (1-nbequ)
#define kde(ee,i) zr(j_kde+(ee-1)*nbequ+i-1)
#define kvi(ee,i) zr(j_kvi+(ee-1)*nbequ+i-1)
!
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_int = sd_int_
    epsi = 100.d0*r8prem()
    errt0 = r8maem()
!
!   1 - Retrieval of the system's state at instant t_i (index=1)
    call intget(sd_int, TIME, iocc=1, rscal=t1, buffer=buffint)
    call intget(sd_int, INDEX, iocc=1, iscal=ind1, buffer=buffint)
    call intget(sd_int, STEP, iocc=1, rscal=dt, buffer=buffint)
    dt2 = dt
!
    call intget(sd_int, DEPL, iocc=1, vr=depl1, lonvec=nbequ, &
                buffer=buffint)
    call intget(sd_int, VITE, iocc=1, vr=vite1, buffer=buffint)
    call intget(sd_int, ACCE, iocc=1, vr=acce1, buffer=buffint)
    call intget(sd_int, FORCE_EX, iocc=1, vr=fext1, buffer=buffint)
!
!   2 - Detection of the initial call to the Runge-Kutta 5/4 algorithm
!       DEPL/2 does not exist in the buffer
    call intget(sd_int, DEPL, iocc=2, lonvec=iret, buffer=buffint)
    if (iret .eq. 0) then
!
!       2.1 - Algorithm initialization
!
        call intinivec(sd_int, PARAMS, 56, vr=par)
!
        call getvr8('SCHEMA_TEMPS', 'TOLERANCE', iocc=1, scal=tol)
        call getvr8('SCHEMA_TEMPS', 'PAS_MINI', iocc=1, scal=dtmin, nbret=iret1)
        call getvr8('SCHEMA_TEMPS', 'PAS_MAXI', iocc=1, scal=dtmax, nbret=iret2)
        call getvr8('SCHEMA_TEMPS', 'ALPHA', iocc=1, scal=alpha, nbret=iret3)
!
        if (iret1 .ne. 1) dtmin = 1.d-10*dt
        if (iret2 .ne. 1) dtmax = 1.d10*dt
        if (iret3 .ne. 1) alpha = 0.d0
!
!       deltadt gives the ratio between dtmin and dtmax, it is considered as an
!       indicator for whether we should adapt or no the time step
        deltadt = abs(dtmax/dtmin-1.d0)
!
!       Dormand-Prince coefficients
!       Note, these are actually saved in PARAMS, a #define is given in the beginning
!       of this subroutine to relate cdp and adp to the PARAMS vector (pointed to by par)
        cdp(1) = 0.0d0
        cdp(2) = 0.2d0
        cdp(3) = 0.3d0
        cdp(4) = 0.8d0
        cdp(5) = 8.d0/9.d0
        cdp(6) = 1.d0
        cdp(7) = 1.d0
!
        adp(1, 1) = 0.2d0
        adp(2, 1) = 3.d0/40.d0
        adp(2, 2) = 9.d0/40.d0
!
        adp(3, 1) = 44.d0/45.d0
        adp(3, 2) = -56.d0/15.d0
        adp(3, 3) = 32.d0/9.d0
!
        adp(4, 1) = 19372.d0/6561.d0
        adp(4, 2) = -25360.d0/2187.d0
        adp(4, 3) = 64448.d0/6561.d0
        adp(4, 4) = -212.d0/729.d0
!
        adp(5, 1) = 9017.d0/3168.d0
        adp(5, 2) = -355.d0/33.d0
        adp(5, 3) = 46732.d0/5247.d0
        adp(5, 4) = 49.d0/176.d0
        adp(5, 5) = -5103.d0/18656.d0
!
        adp(6, 1) = 35.d0/384.d0
        adp(6, 2) = 0.d0
        adp(6, 3) = 500.d0/1113.d0
        adp(6, 4) = 125.d0/192.d0
        adp(6, 5) = -2187.d0/6784.d0
        adp(6, 6) = 11.d0/84.d0
!
        b4b5(1) = 71.0d0/57600.0d0
        b4b5(2) = 0.0d0
        b4b5(3) = -71.0d0/16695.0d0
        b4b5(4) = 71.0d0/1920.0d0
        b4b5(5) = -17253.0d0/339200.0d0
        b4b5(6) = 22.0d0/525.0d0
        b4b5(7) = -1.0d0/40.0d0
!
!       --- Allocate work vectors kde and kdv
        call intinivec(sd_int, WORK1, 6*nbequ, address=j_kde)
        call intinivec(sd_int, WORK2, 6*nbequ, address=j_kvi)
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
        nullify (buffint)
        call intbuff(sd_int, buffint, level=2)
    else
!       --- Algorithm is already initialized, just retrieve DEPL/VITE/ACCE/2
!           Note that these are made identical to index 1 (dcopy before exit)
        call intget(sd_int, DEPL, iocc=2, vr=depl2, buffer=buffint)
        call intget(sd_int, VITE, iocc=2, vr=vite2, buffer=buffint)
        call intget(sd_int, ACCE, iocc=2, vr=acce2, buffer=buffint)
        call intget(sd_int, FORCE_EX, iocc=2, vr=fext2, buffer=buffint)
!
!       --- Retrieve work vectors kde and kdv
        call intget(sd_int, WORK1, address=j_kde, buffer=buffint)
        call intget(sd_int, WORK2, address=j_kvi, buffer=buffint)
!
!       --- Retrieve the algorithm parameters
        call intget(sd_int, PARAMS, vr=par, buffer=buffint)
!
!       --- Retrieve choc parameters save container
        if (nbsavnl .gt. 0) call intget(sd_int, WORK3, vr=nlsav0, buffer=buffint)
!
    end if
!
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
10  continue
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, depl1, b_incx, depl2, b_incy)
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vite1, b_incx, vite2, b_incy)
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, acce1, b_incx, acce2, b_incy)
    if (nbsavnl .gt. 0) then
        b_n = to_blas_int(nbsavnl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, nlsav0, b_incx, nlsav1, b_incy)
    end if
!
!   3 - Loop over all levels (1-6)
    do ee = 1, 6
!
        t2 = t1+dt*cdp(ee+1)
        dt2 = dt*(cdp(ee+1)-cdp(ee))
!
        call intsav(sd_int, TIME, 1, iocc=2, rscal=t2, &
                    buffer=buffint)
        call intsav(sd_int, STEP, 1, iocc=2, rscal=dt2, &
                    buffer=buffint)
!
        b_n = to_blas_int(nbequ)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vite2, b_incx, kde(ee, 1), b_incy)
        b_n = to_blas_int(nbequ)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, acce2, b_incx, kvi(ee, 1), b_incy)
!
        do i = 1, nbequ
            depl2(i) = depl1(i)
            vite2(i) = vite1(i)
            do ss = 1, ee
                depl2(i) = depl2(i)+dt*adp(ee, ss)*kde(ss, i)
                vite2(i) = vite2(i)+dt*adp(ee, ss)*kvi(ss, i)
            end do
        end do
!
        call dtmacce(sd_dtm, sd_int, 2, buffdtm, buffint)
!
    end do
    dt2 = dt
!
!   4 - If an adaptative scheme is requested (dtmax != dtmin)
    if (deltadt .ge. epsi) then
!
!   4.1 - Estimation of the global errors in displacement and velocity
        errd = 0.d0
        errv = 0.d0
        do i = 1, nbequ
!       --- Since b4b5(2) = 0, skip index 2 from the sum-up loop over all 6 levels
            errde = b4b5(1)*kde(1, i)
            errvi = b4b5(1)*kvi(1, i)
            do ee = 3, 6
                errde = errde+b4b5(ee)*kde(ee, i)
                errvi = errvi+b4b5(ee)*kvi(ee, i)
            end do
!
!       --- The last estimation is not saved in kde and kvi, but is present in VITE/ACCE/2
            errde = errde+b4b5(7)*vite2(i)
            errvi = errvi+b4b5(7)*acce2(i)
!
            ddep = abs(depl2(i)-depl1(i))
            dvit = abs(vite2(i)-vite1(i))
!       --- Use the relative tolerance for the local error
            skd = max(tol*abs(depl1(i)+alpha), tol*(ddep+alpha), epsi)
            errd = errd+(errde/skd)**2
!       vitesses. Ajouter atol pour tenir compte d'une tolérance relative.
            skv = max(tol*abs(vite1(i)+alpha), tol*(dvit+alpha), epsi)
            errv = errv+(errvi/skv)**2
        end do
!
!   4.2 - Calculation of a total error, and timestep adaptation
        errt = max(dt*sqrt((errd+errv)/(2*nbequ)), epsi)
!
        seuil1 = (0.9d0/5.0d0)**(6.d0)
        seuil2 = (0.9d0/0.2d0)**(6.d0)
!
!       --- Refining the time step is not helping, the system is probably badly
!           initialized (not in dynamic equilibrium) and the integration can proceed
!           as is toward an equilibrium state
        if ((errt*1.01d0) .gt. errt0) then
            errt = 0.99d0
        end if
!
        if (errt .lt. seuil1) then
            coeff = 5.0d0
        else if (errt .gt. seuil2) then
            coeff = 0.2d0
        else
            coeff = 0.9d0*(errt)**(-1.d0/6.d0)
        end if
        dt2 = min(dtmax, dt*coeff)
!
        if ((errt .ge. 1.d0) .and. (dt2 .ge. dtmin)) then
            dt = dt2
            errt0 = errt
            goto 10
        else if (dt2 .lt. dtmin) then
            call utmess('F', 'ALGORITH5_23')
        end if
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
    call dcopy(b_n, vite2, b_incx, vite1, b_incy)
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, acce2, b_incx, acce1, b_incy)
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, fext2, b_incx, fext1, b_incy)
    call intsav(sd_int, STEP, 1, iocc=1, rscal=dt2, &
                buffer=buffint)
    call intsav(sd_int, TIME, 1, iocc=1, rscal=t1+dt, &
                buffer=buffint)
    call intsav(sd_int, INDEX, 1, iocc=1, iscal=ind1+1, &
                buffer=buffint)
!
    call intsav(sd_int, STEP, 1, iocc=2, rscal=dt, &
                buffer=buffint)
    call intsav(sd_int, TIME, 1, iocc=2, rscal=t1+dt, &
                buffer=buffint)
    call intsav(sd_int, INDEX, 1, iocc=2, iscal=ind1+1, &
                buffer=buffint)
!
!   6 - Set the archiving index to 2
    call intsav(sd_int, IND_ARCH, 1, iscal=2, buffer=buffint)
!
end subroutine
