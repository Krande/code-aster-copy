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
subroutine fft(s, n, ifft)
    implicit none
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/veri32.h"
#include "blas/zscal.h"
!-----------------------------------------------------------------------
! IN,OUT : S    FONCTION A TRANSFORMER
! IN     : N    NOMBRE DE POINTS DE LA FONCTION
! IN     : IFFT > 0 => FFT
!               < 0 => FFT INVERSE
!
    complex(kind=8) :: u, w, t
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
! parametres
    integer(kind=8) :: n, ifft
    complex(kind=8) :: s(n)
! variables locales
    integer(kind=8) :: i, ip, isgn, j, k, l
    integer(kind=8) :: le, le1, m, n2, nm1, nv2
    complex(kind=8) :: calpha
    real(kind=8) :: pi
    blas_int :: b_incx, b_n
!-----------------------------------------------------------------------
!
    m = int(log(dble(n))/log(2.d0))
!    if (m .gt. 30) call veri32()
    n2 = 2**m
    if (n2 .ne. n) then
        m = m+1
!        if (m .gt. 30) call veri32()
        n2 = 2**m
        if (n2 .ne. n) then
            m = m-2
        end if
    end if
    isgn = 1
    if (ifft .lt. 0) isgn = -1
    pi = r8pi()*isgn
    nm1 = n-1
    j = 1
    nv2 = n/2
    do i = 1, nm1
        if (i .ge. j) goto 5
        t = s(j)
        s(j) = s(i)
        s(i) = t
5       continue
        k = nv2
6       continue
        if (k .ge. j) goto 7
        j = j-k
        k = k/2
        goto 6
7       continue
        j = j+k
    end do
    do l = 1, m
!        if (l .gt. 30) call veri32()
        le = 2**l
        le1 = le/2
        u = (1.d0, 0.d0)
        w = dcmplx(cos(-pi/dble(le1)), sin(-pi/dble(le1)))
        do j = 1, le1
            do i = j, n, le
                ip = i+le1
                t = s(ip)*u
                s(ip) = s(i)-t
                s(i) = s(i)+t
            end do
            u = u*w
        end do
    end do
    if (ifft .lt. 0) then
        calpha = dcmplx(1.d0, 0.d0)/(n2*1.d0)
        b_n = to_blas_int(n2)
        b_incx = to_blas_int(1)
        call zscal(b_n, calpha, s, b_incx)
    end if
end subroutine
