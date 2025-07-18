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
subroutine tridia(n, a, lda, d, e, &
                  tau, w)
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/zadder.h"
#include "asterfort/zader2.h"
#include "asterfort/zmvpy.h"
#include "blas/zdotc.h"
    integer(kind=8) :: n, lda
    real(kind=8) :: d(*), e(*)
    complex(kind=8) :: a(lda, *), tau(*), w(*)
!      REDUCTION D'UNE MATRICE HERMITIENNE EN UNE MATRICE TRIDIAGONALE
!                   SYMETRIQUE (METHODE DE HOUSEHOLDER).
!-----------------------------------------------------------------------
! IN  : N    : DIMENSION DE LA MATRICE.
!     : A    : MATRICE HERMITIENNE.
!     : LDA  : DIMENSION DE A.
! OUT : D    : VECTEUR DE DIMENSION N CONTENANT LA DIAGONALE DE LA
!              MATRICE TRIDIAGONALE
!     : E    : VECTEUR DE DIMENSION N CONTENANT LA DIAGONALE
!              SUPERIEURE DE LA MATRICE TRIDIAGONALE DANS E(2:N),
!              E(1) = 0.
!     : TAU  : VECTEUR COMPLEXE DE DIMENSION N CONTENANT LA DIAGONALE
!              DE LA MATRICE UNITAIRE T
!     : W    : VECTEUR COMPLEXE DE DIMENSION N (VECTEUR DE TRAVAIL)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    real(kind=8) :: bb, delta, dgamma, ratio, rho, root, tol, vr
    complex(kind=8) :: temp1
    complex(kind=8) :: vc
    blas_int :: b_incx, b_incy, b_n
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    tol = 0.0d0
    do i = 1, n
        do j = 1, i
            tol = max(abs(dble(a(i, j))), abs(dimag(a(i, j))), tol)
        end do
    end do
    dgamma = r8prem()**2
!
!  --- REALISATION DE N-2 TRANSFORMATIONS SIMILAIRES ---
    do k = 2, n-1
        tau(k) = 0.0d0
        b_n = to_blas_int(n-k+1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        vc = zdotc(b_n, a(k, k-1), b_incx, a(k, k-1), b_incy)
        vr = max(abs(dble(vc)), abs(dimag(vc)))
        if (vr .le. dgamma*tol**2) goto 30
        if (dble(a(k, k-1)) .eq. 0.0d0 .and. dimag(a(k, k-1)) .eq. 0.0d0) then
            a(k, k-1) = sqrt(vr)
            delta = vr
            tau(1) = -a(k, k-1)
        else
            root = abs(a(k, k-1))*sqrt(vr)
            delta = vr+root
            ratio = vr/root
            tau(1) = -ratio*dconjg(a(k, k-1))
            a(k, k-1) = (ratio+1.0d0)*a(k, k-1)
        end if
!
!   --- TRANSFORMATIONS ---
        do j = k, n
            tau(j) = a(j, k-1)/delta
        end do
!
        call zmvpy('LOWER', n-k+1, (1.0d0, 0.0d0), a(k, k), lda, &
                   a(k, k-1), 1, (0.0d0, 0.0d0), w(k), 1)
!                                  RHO = U*NV
        b_n = to_blas_int(n-k+1)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        temp1 = zdotc(b_n, w(k), b_incx, tau(k), b_incy)
        rho = dble(temp1)
        call zader2('LOWER', n-k+1, (-1.0d0, 0.0d0), tau(k), 1, &
                    w(k), 1, a(k, k), lda)
        call zadder('LOWER', n-k+1, rho*delta, tau(k), 1, &
                    a(k, k), lda)
        tau(k) = tau(1)
30      continue
    end do
!
!  --- LA MATRICE A ETE REDUITE EN UNE MATRICE HERMITIENNE
!      TRIDIAGONALE. LA DIAGONALE SUPERIEURE EST TEMPORAIREMENT
!      STOCKEE DANS LE VECTEUR TAU. LA DIAGONALE EST STOCKEE DANS D
    do i = 1, n
        d(i) = dble(a(i, i))
    end do
!
    tau(1) = 1.0d0
    if (n .gt. 1) tau(n) = dconjg(a(n, n-1))
    e(1) = 0.0d0
!
    do i = 2, n
        bb = abs(tau(i))
        e(i) = bb
        a(i, i) = dcmplx(dble(a(i, i)), bb)
        if (bb .eq. 0.0d0) then
            tau(i) = 1.0d0
            bb = 1.0d0
        end if
        tau(i) = tau(i)*tau(i-1)/bb
    end do
!
end subroutine
