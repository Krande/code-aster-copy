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
subroutine zader2(uplo, n, alpha, x, incx, &
                  y, incy, a, lda)
    implicit none
#include "asterf_types.h"
#include "blas/zaxpy.h"
    integer(kind=8) :: n, incx, incy, lda
    complex(kind=8) :: alpha, x(*), y(*), a(lda, *)
    character(len=*) :: uplo
!   CALCUL DE A: MATRICE HERMITIENNE
!   A = A + ALPHA*X*CONJG(Y)' + CONJG(ALPHA)*Y*CONJG(X)'
!-----------------------------------------------------------------------
! IN  : UPLO : INDIQUE LE MODE DE STOCKAGE DE LA MATRICE.
!              SI UPLO EST 'U' ALORS SEULEMENT LA PARTIE SUPERIEURE DE A
!              EST UTILISEE. SI UPLO EST 'L', ALORS LA PARTIE INFERIEURE
!              EST UTILISEE.
!     : N    : DIMENSION DE LA MATRICE A.
!     : ALPHA: SCALAIRE.
!     : X    : DVECTEURE COMPLEXE DE LONGUEUR (N-1)*IABS(INCX)+1.
!     : INCX : DEPLACEMENT ENTRE LES ELEMENTS DE X.
!     : Y    : DVECTEURE COMPLEXE DE LONGUEUR (N-1)*IABS(INCY)+1.
!     : INCY : DEPLACEMENT ENTRE LES ELEMENTS DE Y.
! I/O : A    : MATRICE COMPLEXE DE DIMENSION N.
! IN  : LDA  : DIMENSION DE A
!-----------------------------------------------------------------------
    integer(kind=8) :: ix, iy, j
    complex(kind=8) :: tempx, tempy, temp1
    aster_logical :: upper
    real(kind=8) :: dble
    blas_int :: b_incx, b_incy, b_n
!
    if (n .eq. 0 .or. alpha .eq. (0.0d0, 0.0d0)) goto 999
!
    ix = 1
    iy = 1
    if (incx .lt. 0) ix = 1-(n-1)*incx
    if (incy .lt. 0) iy = 1-(n-1)*incy
!
    upper = (uplo(1:1) .eq. 'U') .or. (uplo(1:1) .eq. 'u')
!
    do j = 1, n
        tempx = dconjg(alpha*x(ix))
        tempy = alpha*dconjg(y(iy))
        if (upper) then
            if (incx .ge. 0) then
                b_n = to_blas_int(j-1)
                b_incx = to_blas_int(incx)
                b_incy = to_blas_int(1)
                call zaxpy(b_n, tempy, x, b_incx, a(1, j), &
                           b_incy)
            else
                b_n = to_blas_int(j-1)
                b_incx = to_blas_int(incx)
                b_incy = to_blas_int(1)
                call zaxpy(b_n, tempy, x(ix-incx), b_incx, a(1, j), &
                           b_incy)
            end if
            if (incy .ge. 0) then
                b_n = to_blas_int(j-1)
                b_incx = to_blas_int(incy)
                b_incy = to_blas_int(1)
                call zaxpy(b_n, tempx, y, b_incx, a(1, j), &
                           b_incy)
            else
                b_n = to_blas_int(j-1)
                b_incx = to_blas_int(incy)
                b_incy = to_blas_int(1)
                call zaxpy(b_n, tempx, y(iy-incy), b_incx, a(1, j), &
                           b_incy)
            end if
        else
            if (incx .ge. 0) then
                b_n = to_blas_int(n-j)
                b_incx = to_blas_int(incx)
                b_incy = to_blas_int(1)
                call zaxpy(b_n, tempy, x(ix+incx), b_incx, a(j+1, j), &
                           b_incy)
            else
                b_n = to_blas_int(n-j)
                b_incx = to_blas_int(incx)
                b_incy = to_blas_int(1)
                call zaxpy(b_n, tempy, x, b_incx, a(j+1, j), &
                           b_incy)
            end if
            if (incy .ge. 0) then
                b_n = to_blas_int(n-j)
                b_incx = to_blas_int(incy)
                b_incy = to_blas_int(1)
                call zaxpy(b_n, tempx, y(iy+incy), b_incx, a(j+1, j), &
                           b_incy)
            else
                b_n = to_blas_int(n-j)
                b_incx = to_blas_int(incy)
                b_incy = to_blas_int(1)
                call zaxpy(b_n, tempx, y, b_incx, a(j+1, j), &
                           b_incy)
            end if
        end if
        temp1 = a(j, j)+y(iy)*tempx+x(ix)*tempy
        a(j, j) = dble(temp1)
        ix = ix+incx
        iy = iy+incy
    end do
!
999 continue
end subroutine
