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
subroutine intnewm(sd_dtm_, sd_int_, buffdtm, buffint)
    implicit none
!
! person_in_charge: hassan.berro at edf.fr
!
! intnewm : Integrate from t_i to t_i+1 the differential equations of motion
!           using the NEWMARK integration method.
!
#include "jeveux.h"
#include "blas/dcopy.h"
#include "asterc/r8prem.h"
#include "asterfort/dtmforc.h"
#include "asterfort/getvr8.h"
#include "asterfort/intbuff.h"
#include "asterfort/intget.h"
#include "asterfort/intnewm_oper.h"
#include "asterfort/intinivec.h"
#include "asterfort/intsav.h"
#include "asterfort/pmavec.h"
#include "asterfort/rrlds.h"
!
!   -0.1- Input/output arguments
    character(len=*), intent(in) :: sd_dtm_
    character(len=*), intent(in) :: sd_int_
    integer(kind=8), pointer :: buffdtm(:)
    integer(kind=8), pointer :: buffint(:)
!
!   -0.2- Local variables
    integer(kind=8) :: i, nbequ, ind1, iret
    integer(kind=8) :: upmat
    real(kind=8) :: t1, dt, bet, gam, coeff, epsi
    real(kind=8) :: dtold
    character(len=8) :: sd_dtm, sd_int
!
    real(kind=8), pointer :: depl1(:) => null()
    real(kind=8), pointer :: vite1(:) => null()
    real(kind=8), pointer :: acce1(:) => null()
    real(kind=8), pointer :: fext1(:) => null()
    real(kind=8), pointer :: depl2(:) => null()
    real(kind=8), pointer :: vite2(:) => null()
    real(kind=8), pointer :: acce2(:) => null()
!
    real(kind=8), pointer :: mgen(:) => null()
    real(kind=8), pointer :: kgen(:) => null()
    real(kind=8), pointer :: agen(:) => null()
!
    real(kind=8), pointer :: par(:) => null()
    real(kind=8), pointer :: ktilda(:) => null()
    real(kind=8), pointer :: ftild1(:) => null()
    real(kind=8), pointer :: ftild2(:) => null()
    real(kind=8), pointer :: ftild3(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
#define a(k) par(k)
#define norm_coef par(9)
#define beta par(10)
#define gamma par(11)
#define mdiag_r par(12)
#define kdiag_r par(13)
#define cdiag_r par(14)
!
#define mdiag (nint(mdiag_r).eq.1)
#define kdiag (nint(kdiag_r).eq.1)
#define cdiag (nint(cdiag_r).eq.1)
!
#define m(row,col) mgen((col-1)*nbequ+row)
#define k(row,col) kgen((col-1)*nbequ+row)
#define c(row,col) agen((col-1)*nbequ+row)
!
#define kt(row,col) ktilda((col-1)*nbequ+row)
#define ft1(row,col) ftild1((col-1)*nbequ+row)
#define ft2(row,col) ftild2((col-1)*nbequ+row)
#define ft3(row,col) ftild3((col-1)*nbequ+row)
!
!
!   0 - Initializations
    sd_dtm = sd_dtm_
    sd_int = sd_int_
    epsi = r8prem()
    coeff = 1.d0
!
!   1 - Retrieval of the system's state at instant t_i (index=1)
    call intget(sd_int, TIME, iocc=1, rscal=t1, buffer=buffint)
    call intget(sd_int, INDEX, iocc=1, iscal=ind1, buffer=buffint)
    call intget(sd_int, STEP, iocc=1, rscal=dt, buffer=buffint)
!
    call intget(sd_int, DEPL, iocc=1, vr=depl1, lonvec=nbequ, &
                buffer=buffint)
    call intget(sd_int, VITE, iocc=1, vr=vite1, buffer=buffint)
    call intget(sd_int, ACCE, iocc=1, vr=acce1, buffer=buffint)
!
!   2 - Detection of the initial call to the Newmark algorithm
!       DEPL/2 does not exist in the buffer
    call intget(sd_int, DEPL, iocc=2, lonvec=iret, buffer=buffint)
    if (iret .eq. 0) then
!
!       2.1 - Algorithm initialization :
!                   - Retrieval and saving the scheme's parameters : beta and gamma
!                   - Calculating all necessary coefficients and operators
!
        call getvr8('SCHEMA_TEMPS', 'BETA', iocc=1, scal=bet)
        call getvr8('SCHEMA_TEMPS', 'GAMMA', iocc=1, scal=gam)
!
!       --- Algorithm parameters :
!       [1-8]
        call intinivec(sd_int, PARAMS, 14, vr=par)
!
        beta = bet
        gamma = gam
!
        a(1) = 1.d0/(beta*dt*dt)
        a(2) = gamma/(beta*dt)
        a(3) = 1.d0/(beta*dt)
        a(4) = (1.d0/(2*beta))-1
        a(5) = (gamma/beta)-1
        a(6) = (dt/2)*((gamma/beta)-2)
        a(7) = dt*(1-gamma)
        a(8) = gamma*dt
!
!       --- Normalizing coefficient for ktilda, ftilda*, for reducing numerical errors
!           encountered with tiny time steps
        norm_coef = 1.d0
!
!       --- Detection of the matrix types (full/diagonal) for M, K, and C
        mdiag_r = 0.d0
        call intget(sd_int, MASS_FUL, iocc=1, lonvec=iret, buffer=buffint)
        if (iret .gt. 0) then
            call intget(sd_int, MASS_FUL, iocc=1, vr=mgen, buffer=buffint)
        else
            call intget(sd_int, MASS_DIA, iocc=1, vr=mgen, buffer=buffint)
            mdiag_r = 1.d0
        end if
!
        kdiag_r = 0.d0
        call intget(sd_int, RIGI_FUL, iocc=1, lonvec=iret, buffer=buffint)
        if (iret .gt. 0) then
            call intget(sd_int, RIGI_FUL, iocc=1, vr=kgen, buffer=buffint)
        else
            call intget(sd_int, RIGI_DIA, iocc=1, vr=kgen, buffer=buffint)
            kdiag_r = 1.d0
        end if
!
        cdiag_r = 0.d0
        call intget(sd_int, AMOR_FUL, iocc=1, lonvec=iret, buffer=buffint)
        if (iret .gt. 0) then
            call intget(sd_int, AMOR_FUL, iocc=1, vr=agen, buffer=buffint)
        else
            call intget(sd_int, AMOR_DIA, iocc=1, vr=agen, buffer=buffint)
            cdiag_r = 1.d0
        end if
!
!       --- Memory allocation for  operators ktilda, ftild1, ftild2, and ftild3
        if (mdiag .and. kdiag .and. cdiag) then
!           --- Operator ktilda is diagonal
            call intinivec(sd_int, WORK1, nbequ, vr=ktilda)
        else
!           --- Operator ktilda is full
            call intinivec(sd_int, WORK1, nbequ*nbequ, vr=ktilda)
        end if
!
        if (mdiag .and. cdiag) then
!           --- Operators ftild[1,2,3] are diagonal
            call intinivec(sd_int, WORK2, nbequ, vr=ftild1)
            call intinivec(sd_int, WORK3, nbequ, vr=ftild2)
            call intinivec(sd_int, WORK4, nbequ, vr=ftild3)
        else
!           --- Operators ftild[1,2,3] are then full
            call intinivec(sd_int, WORK2, nbequ*nbequ, vr=ftild1)
            call intinivec(sd_int, WORK3, nbequ*nbequ, vr=ftild2)
            call intinivec(sd_int, WORK4, nbequ*nbequ, vr=ftild3)
        end if
!
!       --- Calculation of ktilda, ftild1, ftild2, and ftild3 according to :
!           kt (i,j) = a1*m(i,j) + k(i,j) + a2*c(i,j)
!           ft1(i,j) = a3*m(i,j) + a5*c(i,j)
!           ft2(i,j) = a1*m(i,j) + a2*c(i,j)
!           ft3(i,j) = a4*m(i,j) + a6*c(i,j)
        call intnewm_oper(nbequ, par, mgen, kgen, agen, &
                          ktilda, ftild1, ftild2, ftild3)
!
!       --- Allocate DEPL/VITE/ACCE/2 (t_i+1)
        call intinivec(sd_int, DEPL, nbequ, iocc=2, vr=depl2)
        call intinivec(sd_int, VITE, nbequ, iocc=2, vr=vite2)
        call intinivec(sd_int, ACCE, nbequ, iocc=2, vr=acce2)
!
        dtold = dt
!
        nullify (buffint)
        call intbuff(sd_int, buffint, level=2)
!
    else
!       --- Algorithm is already initialized, retrieval of all operators from work vectors
        call intget(sd_int, WORK1, vr=ktilda, buffer=buffint)
        call intget(sd_int, WORK2, vr=ftild1, buffer=buffint)
        call intget(sd_int, WORK3, vr=ftild2, buffer=buffint)
        call intget(sd_int, WORK4, vr=ftild3, buffer=buffint)
!
        call intget(sd_int, PARAMS, vr=par, buffer=buffint)
!
!       --- Retrieval of already allocated DEPL/VITE/ACCE/2 (t_i+1)
        call intget(sd_int, DEPL, iocc=2, vr=depl2, buffer=buffint)
        call intget(sd_int, VITE, iocc=2, vr=vite2, buffer=buffint)
        call intget(sd_int, ACCE, iocc=2, vr=acce2, buffer=buffint)
!
        call intget(sd_int, STEP, iocc=2, rscal=dtold, buffer=buffint)
    end if
!
!   3 - Computing FORCE/2 (t_i+1)
    call intsav(sd_int, TIME, 1, iocc=2, rscal=t1+dt, &
                buffer=buffint)
    call intsav(sd_int, STEP, 1, iocc=2, rscal=dt, &
                buffer=buffint)
    call intsav(sd_int, INDEX, 1, iocc=2, iscal=ind1+1, &
                buffer=buffint)
    call dtmforc(sd_dtm, sd_int, 2, buffdtm, buffint)
    call intget(sd_int, FORCE_EX, iocc=2, vr=fext1, buffer=buffint)
!
    coeff = dt/dtold
!
    call intget(sd_int, MAT_UPDT, iscal=upmat, buffer=buffint)
!
!   4 - Updating the operators in the event of a change in dt or a matrix update
!       (case of implicit consideration of non-linearities)
    if ((upmat .eq. 1) .or. (abs(coeff-1.d0) .ge. epsi)) then
        if (upmat .eq. 1) then
            mdiag_r = 0.d0
            call intget(sd_int, MASS_FUL, iocc=1, lonvec=iret, buffer=buffint)
            if (iret .gt. 0) then
                call intget(sd_int, MASS_FUL, iocc=1, vr=mgen, buffer=buffint)
            else
                call intget(sd_int, MASS_DIA, iocc=1, vr=mgen, buffer=buffint)
                mdiag_r = 1.d0
            end if
            !
            kdiag_r = 0.d0
            call intget(sd_int, RIGI_FUL, iocc=1, lonvec=iret, buffer=buffint)
            if (iret .gt. 0) then
                call intget(sd_int, RIGI_FUL, iocc=1, vr=kgen, buffer=buffint)
            else
                call intget(sd_int, RIGI_DIA, iocc=1, vr=kgen, buffer=buffint)
                kdiag_r = 1.d0
            end if
            !
            cdiag_r = 0.d0
            call intget(sd_int, AMOR_FUL, iocc=1, lonvec=iret, buffer=buffint)
            if (iret .gt. 0) then
                call intget(sd_int, AMOR_FUL, iocc=1, vr=agen, buffer=buffint)
            else
                call intget(sd_int, AMOR_DIA, iocc=1, vr=agen, buffer=buffint)
                cdiag_r = 1.d0
            end if
!
            call intsav(sd_int, MAT_UPDT, 1, iscal=0, buffer=buffint)
        end if
!
!       Updating the Newmark parameters and integration operators
        a(1) = 1.d0/(beta*dt*dt)
        a(2) = gamma/(beta*dt)
        a(3) = 1.d0/(beta*dt)
        a(4) = (1.d0/(2*beta))-1
        a(5) = (gamma/beta)-1
        a(6) = (dt/2)*((gamma/beta)-2)
        a(7) = dt*(1-gamma)
        a(8) = gamma*dt
!
        call intnewm_oper(nbequ, par, mgen, kgen, agen, &
                          ktilda, ftild1, ftild2, ftild3)
    end if
!
!   5 - Calculating DEPL/2 (t_i+1)according to the equation
!      ktilda(i,j) * depl2(i) =  fext1(i)              +  ftild1(i,j)*vite1(i)  + &
!                                ftild2(i,j)*depl1(i)  +  ftild3(i,j)*acce1(i)
!
!   --- Right hand side vector calculation
    if (size(ftild1) .eq. nbequ) then
        do i = 1, nbequ
            fext1(i) = fext1(i)/norm_coef+ftild1(i)*vite1(i)+ftild2(i)*depl1(i)+ftild3(i)*acce1(&
                       &i)
        end do
    else
        do i = 1, nbequ
            fext1(i) = fext1(i)/norm_coef
        end do
        call pmavec('CUMUL', nbequ, ftild1, vite1, fext1)
        call pmavec('CUMUL', nbequ, ftild2, depl1, fext1)
        call pmavec('CUMUL', nbequ, ftild3, acce1, fext1)
    end if
    b_n = to_blas_int(nbequ)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, fext1, b_incx, depl2, b_incy)
!
!   --- Linear equation AX=B resolution
    if (size(ktilda) .eq. nbequ) then
        do i = 1, nbequ
            depl2(i) = depl2(i)/ktilda(i)
        end do
    else
        call rrlds(ktilda, nbequ, nbequ, depl2, 1)
    end if
!
!   6 - Calculating VITE/ACCE/2 (t_i+1) according to the equation
!       acce2(i) = -a4*acce1(i) + a1*(depl2(i) - depl1(i) - dt*vite1(i))
!       vite2(i) =     vite1(i) + a7*acce1(i) + a8*acce2(i)
!
    do i = 1, nbequ
        acce2(i) = -a(4)*acce1(i)+a(1)*(depl2(i)-depl1(i)-dt*vite1(i))
        vite2(i) = vite1(i)+a(7)*acce1(i)+a(8)*acce2(i)
    end do
!
!   7 - Preparing the algorithm for the next step, copy index 2 in 1
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
    call intsav(sd_int, STEP, 1, iocc=1, rscal=dt, &
                buffer=buffint)
    call intsav(sd_int, TIME, 1, iocc=1, rscal=t1+dt, &
                buffer=buffint)
    call intsav(sd_int, INDEX, 1, iocc=1, iscal=ind1+1, &
                buffer=buffint)
!
!   8 - Set the archiving index to 2
    call intsav(sd_int, IND_ARCH, 1, iscal=2, buffer=buffint)
!
end subroutine
