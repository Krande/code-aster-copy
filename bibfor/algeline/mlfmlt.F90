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
subroutine mlfmlt(b, f, y, ldb, n, &
                  p, l, opta, optb, nb)
!
!     B = B - F*Y PAR BLOCS
    implicit none
#include "blas/dgemm.h"
    integer(kind=8) :: ldb, n, p, l, nb
    real(kind=8) :: b(ldb, l), f(n, p), y(ldb, l)
    integer(kind=8) :: m, npb, restp, nlb, restl
    integer(kind=8) :: i, j, ib, jb
    real(kind=8) :: beta, alpha
    integer(kind=8) :: opta, optb
    character(len=1) :: tra, trb
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
!
    tra = 'T'
    if (opta .eq. 1) tra = 'N'
    trb = 'T'
    if (optb .eq. 1) trb = 'N'
    alpha = -1.d0
    beta = 1.d0
    m = n-p
    npb = p/nb
    nlb = l/nb
    restp = p-(nb*npb)
    restl = l-(nb*nlb)
    if (nlb .gt. 0) then
        do j = 1, nlb
            jb = nb*(j-1)+1
            do i = 1, npb
                ib = nb*(i-1)+1
                b_ldc = to_blas_int(ldb)
                b_ldb = to_blas_int(ldb)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(nb)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(m)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, f(1, ib), b_lda, y(1, jb), b_ldb, &
                           beta, b(ib, jb), b_ldc)
            end do
            if (restp .gt. 0) then
                ib = nb*npb+1
                b_ldc = to_blas_int(ldb)
                b_ldb = to_blas_int(ldb)
                b_lda = to_blas_int(n)
                b_m = to_blas_int(restp)
                b_n = to_blas_int(nb)
                b_k = to_blas_int(m)
                call dgemm(tra, trb, b_m, b_n, b_k, &
                           alpha, f(1, ib), b_lda, y(1, jb), b_ldb, &
                           beta, b(ib, jb), b_ldc)
            end if
        end do
    end if
    if (restl .gt. 0) then
        jb = nb*nlb+1
        do i = 1, npb
            ib = nb*(i-1)+1
            b_ldc = to_blas_int(ldb)
            b_ldb = to_blas_int(ldb)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(nb)
            b_n = to_blas_int(restl)
            b_k = to_blas_int(m)
            call dgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, f(1, ib), b_lda, y(1, jb), b_ldb, &
                       beta, b(ib, jb), b_ldc)
        end do
        if (restp .gt. 0) then
            ib = nb*npb+1
            b_ldc = to_blas_int(ldb)
            b_ldb = to_blas_int(ldb)
            b_lda = to_blas_int(n)
            b_m = to_blas_int(restp)
            b_n = to_blas_int(restl)
            b_k = to_blas_int(m)
            call dgemm(tra, trb, b_m, b_n, b_k, &
                       alpha, f(1, ib), b_lda, y(1, jb), b_ldb, &
                       beta, b(ib, jb), b_ldc)
        end if
    end if
end subroutine
