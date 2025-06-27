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
subroutine zmulmv(trans, m, n, alpha, a, &
                  lda, x, incx, beta, y, &
                  incy)
    implicit none
#include "asterfort/vecinc.h"
#include "asterfort/zmult.h"
#include "blas/zaxpy.h"
#include "blas/zdotc.h"
#include "blas/zdotu.h"
    integer(kind=8) :: m, n, lda, incx, incy
    complex(kind=8) :: alpha, beta, x(*), y(*)
    character(len=*) :: trans
!  CALCUL DU PRODUIT D'UNE MATRICE PAR UN VECTEUR (OPTION 'N' 'T' 'C').
!-----------------------------------------------------------------------
! IN  : TRANS: CARACTERE SPECIFIANT L'OPERATION A REALISER.
!                 TRANS               OPERATION
!              'N'             Y = ALPHA*A*X + BETA*Y
!              'T'             Y = ALPHA*A'*X + BETA*Y
!              'C'             Y = ALPHA*CONJG(A)'*X + BETA*Y
!     : M    : NOMBRE DE LIGNES DE A.
!     : N    : NOMBRE DE COLONNES DE A.
!     : ALPHA: SCALAIRE.
!     : A    : MATRICE COMPLEXE DE DIMENSION M*N
!     : LDA  : DIMENSION DE A
!     : X    : VECTEUR COMLEXE DE LONGUEUR (N-1)*IABS(INCX)+1 LORSQUE
!              TRANS EST EGAL A 'N' ET DE LONGUEUR (M-1)*IABS(INCX)+1
!              SINON.
!     : INCX : DEPLACEMENT ENTRE LES ELEMENTS DE X.
!     : BETA : COMPLEXE.'LORSQUE BETA EGAL ZERO, Y EST NON CALCULE.
! I/O : Y    :  (N-1)*IABS(INCY)+1
!               (M-1)*IABS(INCY)+1
! OUT : INCY : DEPLACEMENT ENTRE LES ELEMENTS DE Y.
!-----------------------------------------------------------------------
!                                  SPECIFICATIONS FOR LOCAL VARIABLES
    integer(kind=8) :: i, ix, iy, ky, lenx, leny
    complex(kind=8) :: a(*)
    integer(kind=8) :: kx
    blas_int :: b_incx, b_incy, b_n
!
    if (m .eq. 0 .or. n .eq. 0 .or. alpha .eq. (0.0d0, 0.0d0) .and. beta .eq. (1.0d0, 0.0d0)) &
        goto 999
!
    if (trans(1:1) .eq. 'N' .or. trans(1:1) .eq. 'n') then
        lenx = n
        leny = m
    else
        lenx = m
        leny = n
    end if
!
    ix = 1
    iy = 1
    if (incx .lt. 0) ix = (-lenx+1)*incx+1
    if (incy .lt. 0) iy = (-leny+1)*incy+1
!
    if (beta .eq. (1.0d0, 0.0d0)) then
    else if (incy .eq. 0) then
        if (beta .eq. (0.0d0, 0.0d0)) then
            y(1) = (0.0d0, 0.0d0)
        else
            y(1) = beta**leny*y(1)
        end if
    else if (beta .eq. (0.0d0, 0.0d0)) then
        call vecinc(leny, (0.0d0, 0.0d0), y, inc=abs(incy))
    else
        call zmult(leny, beta, y, abs(incy))
    end if
!
    if (alpha .eq. (0.0d0, 0.0d0)) goto 999
!
    if (trans(1:1) .eq. 'N' .or. trans(1:1) .eq. 'n') then
        kx = ix
        do i = 1, n
            b_n = to_blas_int(m)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(incy)
            call zaxpy(b_n, alpha*x(kx), a(lda*(i-1)+1), b_incx, y, &
                       b_incy)
            kx = kx+incx
        end do
    else if (trans(1:1) .eq. 'T' .or. trans(1:1) .eq. 't') then
!
        ky = iy
        do i = 1, n
            b_n = to_blas_int(m)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(incx)
            y(ky) = y(ky)+alpha*zdotu(b_n, a(lda*(i-1)+1), b_incx, x, b_incy)
            ky = ky+incy
        end do
!
    else
        ky = iy
        do i = 1, n
            b_n = to_blas_int(m)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(incx)
            y(ky) = y(ky)+alpha*zdotc(b_n, a(lda*(i-1)+1), b_incx, x, b_incy)
            ky = ky+incy
        end do
    end if
!
999 continue
end subroutine
