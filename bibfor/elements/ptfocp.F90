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
subroutine ptfocp(itype, option, xl, nno, nc, &
                  pgl, fer, fei)
!
!
! ------------------------------------------------------------------------------
!
!
! ------------------------------------------------------------------------------
!
    implicit none
!
    integer(kind=8) :: itype, nno, nc
    character(len=*) :: option
    real(kind=8) :: fer(12), fei(12), pgl(3, 3)
    real(kind=8) :: xl
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/provec.h"
#include "asterfort/pscvec.h"
#include "asterfort/ptfop1.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "blas/ddot.h"
!
! ------------------------------------------------------------------------------
!
    integer(kind=8) :: i, icoec, icoer, iret, lforc, lx, ncc, nnoc
    real(kind=8) :: coef1, coef2, s, s2, xxx, s3, s4, s5
    real(kind=8) :: u(3), v(3), w(6), w2(3)
    real(kind=8) :: qr(12), qqr(12), qi(12), qqi(12)
    character(len=16) :: ch16
    aster_logical :: global, normal
    blas_int :: b_incx, b_incy, b_n
!
! ------------------------------------------------------------------------------
!
    qr(:) = 0.0d0
    qqr(:) = 0.0d0
    fer(1:12) = 0.0d0
    qi(:) = 0.0d0
    qqi(:) = 0.0d0
    fei(1:12) = 0.0d0
    nnoc = 1
    ncc = 6
    global = .false.
!
    if (option .eq. 'CHAR_MECA_FC1D1D') then
        call jevech('PGEOMER', 'L', lx)
        lx = lx-1
        do i = 1, 3
            w(i) = zr(lx+i)
            w(i+3) = zr(lx+i+3)
            w2(i) = w(i+3)-w(i)
        end do
!
!       force poutre à valeurs complexes
        call jevech('PFC1D1D', 'L', lforc)
        xxx = abs(dble(zc(lforc+6)))
        global = xxx .lt. 1.d-3
        normal = xxx .gt. 1.001d0
        do i = 1, 3
            qr(i) = dble(zc(lforc-1+i))
            qr(i+6) = qr(i)
            qi(i) = dimag(zc(lforc-1+i))
            qi(i+6) = qi(i)
            xxx = abs(dble(zc(lforc-1+3+i)))
            if (xxx .gt. 1.d-20) then
                call utmess('F', 'ELEMENTS2_46')
            end if
        end do
!
        if (normal) then
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, w2, b_incx, w2, b_incy)
            s2 = 1.d0/s
            call provec(w2, qr(1), u)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, u, b_incx, u, b_incy)
            s3 = sqrt(s)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, qr(1), b_incx, qr(1), b_incy)
            s4 = sqrt(s)
            s5 = s3*sqrt(s2)/s4
            call provec(u, w2, v)
            call pscvec(3, s2, v, u)
            call pscvec(3, s5, u, qr(1))
            call provec(w2, qr(7), u)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, u, b_incx, u, b_incy)
            s3 = sqrt(s)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, qr(7), b_incx, qr(7), b_incy)
            s4 = sqrt(s)
            s5 = s3*sqrt(s2)/s4
            call provec(u, w2, v)
            call pscvec(3, s2, v, u)
            call pscvec(3, s5, u, qr(7))
!
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, w2, b_incx, w2, b_incy)
            s2 = 1.d0/s
            call provec(w2, qi(1), u)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, u, b_incx, u, b_incy)
            s3 = sqrt(s)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, qi(1), b_incx, qi(1), b_incy)
            s4 = sqrt(s)
            s5 = s3*sqrt(s2)/s4
            call provec(u, w2, v)
            call pscvec(3, s2, v, u)
            call pscvec(3, s5, u, qi(1))
            call provec(w2, qi(7), u)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, u, b_incx, u, b_incy)
            s3 = sqrt(s)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            s = ddot(b_n, qi(7), b_incx, qi(7), b_incy)
            s4 = sqrt(s)
            s5 = s3*sqrt(s2)/s4
            call provec(u, w2, v)
            call pscvec(3, s2, v, u)
            call pscvec(3, s5, u, qi(7))
        end if
!       passage repere local du vecteur force (si necessaire)
        if (global .or. normal) then
            call utpvgl(nno, nc, pgl, qr(1), qqr(1))
            call utpvgl(nno, nc, pgl, qi(1), qqi(1))
        else
            do i = 1, 12
                qqr(i) = qr(i)
                qqi(i) = qi(i)
            end do
        end if
!       a cause des chargements variables
        coef1 = 1.0d0
        coef2 = 1.0d0
!
    else
        ch16 = option
        call utmess('F', 'ELEMENTS2_47', sk=ch16)
    end if
!
! ------------------------------------------------------------------------------
!   recuperation du coef_mult
    call tecach('NNO', 'PCOEFFR', 'L', iret, iad=icoer)
    call tecach('NNO', 'PCOEFFC', 'L', iret, iad=icoec)
!
    if (icoer .ne. 0) then
        do i = 1, 12
            qqr(i) = qqr(i)*zr(icoer)
            qqi(i) = qqi(i)*zr(icoer)
        end do
    else if (icoec .ne. 0) then
        do i = 1, 12
            qqr(i) = qqr(i)*dble(zc(icoec))
            qqi(i) = qqi(i)*dimag(zc(icoec))
        end do
    end if
!
    call ptfop1(itype, nc, coef1, coef2, xl, &
                qqr, fer)
    call ptfop1(itype, nc, coef1, coef2, xl, &
                qqi, fei)
!
end subroutine
