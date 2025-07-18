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
subroutine diatri(n, d, e, vector, evec, &
                  ldevec)
    implicit none
#include "asterf_types.h"
#include "asterc/r8miem.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/r8rotg.h"
#include "asterfort/r8sqrt.h"
#include "blas/dcopy.h"
#include "blas/drot.h"
#include "blas/dscal.h"
#include "blas/dswap.h"
#include "blas/idamax.h"
    integer(kind=8) :: n, ldevec
    real(kind=8) :: d(*), e(*), evec(ldevec, *)
    aster_logical :: vector
!   CALCUL DES VALEURS PROPRES ET VECTEURS PROPRES (OPTION => "VECTOR")
!      SUITE A LA TRANSFORMATION DE HOUSEHOLDER ( SUB. TRIDIA).
!-----------------------------------------------------------------------
! IN  : N    : DIMENSION DES MATRICES.
! I/O : D    : VECTEUR DE REELS DE LONGUEUR N.
!         IN : CONTIENT LA DIAGONALE DE LA MATRICE.
!        OUT : CONTIENT LES VALEURS PROPRES DANS L'ORDRE CROISSANT.
!     : E    : VECTEUR DE REELS DE LONGUEUR N.
!         IN : CONTIENT LES ELEMENTS DE LA DIAGONALE, E(1) ARBITRAIRE
!        OUT : VALEURS QUELCONQUES.
! IN  :VECTOR: VARIABLE LOGIQUE  .TRUE. SI LES VECTEURS PROPRES SONT
!              CALCULEES.
! I/O : EVEC : MATRICE REELLE D'ORDRE N.
!         IN : MATRICE TRANSFORMEE UTILISEE POUR TRANFORMER LA MATRICE
!              INITIALE EN UNE MATRICE TRIDIAGONALE.SI LES VECTEURS
!              PROPRES SONT CALCULES, ELLE CONTIENT LA MATRICE IDENTITE.
!        OUT : LE VECTEUR PROPRE ASSOCIE A EVAL(J) (J-TH VALEUR PROPRE
!              EST STOCKE DANS LA J-TH CLONNE.
! IN  :LDEVEC: DIMENSION EXACTE DE EVEC.
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iter, j, k, l, m
    real(kind=8) :: b, c, f, g, p, r, s, scale, tiny, tol
    blas_int :: b_incx, b_incy, b_n
!
    if (n .eq. 1) goto 9000
!
    b_n = to_blas_int(n-1)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, e(2), b_incx, e(1), b_incy)
    e(n) = 0.0d0
!
    tiny = 100.0d0*r8miem()
    tol = r8prem()
    iter = 0
    do l = 1, n
!    --- RECHERCHE DE LA PLUS PETITE VALEUR DE LA DIAGONALE SUPERIEURE.
10      continue
        do m = l, n
            if (m .eq. n) goto 30
            if (abs(e(m)) .le. max(tol*(abs(d(m))+abs(d(m+1))), tiny)) goto 30
        end do
!
30      continue
        p = d(l)
        if (m .eq. l) goto 60
        if (iter .eq. 30*n) then
!C            WRITE(6,*)  'THE ITERATION FOR THE EIGENVALUES DID '//
!C     &                  'NOT CONVERGE.'
            ASSERT(.false.)
        end if
        iter = iter+1
!       --- VALEUR DE SHIFT ---
        g = (d(l+1)-p)/(2.0d0*e(l))
        r = r8sqrt(g, 1.0d0)
        g = d(m)-p+e(l)/(g+sign(r, g))
        s = 1.0d0
        c = 1.0d0
        p = 0.0d0
!
        do i = m-1, l, -1
            f = s*e(i)
            b = c*e(i)
            call r8rotg(g, f, c, s)
            e(i+1) = g
            if (g .eq. 0.0d0) goto 50
            g = d(i+1)-p
            r = (d(i)-g)*s+2.0d0*c*b
            p = s*r
            d(i+1) = g+p
            g = c*r-b
!
            if (vector) then
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call drot(b_n, evec(1, i+1), b_incx, evec(1, i), b_incy, &
                          c, s)
            end if
!
        end do
!
        d(l) = d(l)-p
        e(l) = g
        e(m) = 0.0d0
        goto 10
!
50      continue
        d(i+1) = d(i+1)-p
        e(m) = 0.0d0
        goto 10
60      continue
    end do
!    --- POSITION DES VALEURS ET VECTERUS PROPRES ---
    do i = 1, n-1
        k = i
        p = d(i)
!
        do j = i+1, n
            if (d(j) .lt. p) then
                k = j
                p = d(j)
            end if
        end do
!
        if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            if (vector) then
                b_n = to_blas_int(n)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dswap(b_n, evec(1, i), b_incx, evec(1, k), b_incy)
            end if
        end if
!
    end do
!          --- NORMALISATION DES VECTEURS PROPRES ---
    if (vector) then
        do j = 1, n
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            i = idamax(b_n, evec(1, j), b_incx)
            scale = evec(i, j)
            b_n = to_blas_int(n)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.0d0/scale, evec(1, j), b_incx)
        end do
    end if
!
9000 continue
end subroutine
