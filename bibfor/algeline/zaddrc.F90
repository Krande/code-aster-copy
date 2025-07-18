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
subroutine zaddrc(m, n, alpha, x, incx, &
                  y, incy, a, lda)
    implicit none
#include "blas/zaxpy.h"
    integer(kind=8) :: m, n, incx, incy, lda
    complex(kind=8) :: alpha, x(*), y(*)
    complex(kind=8) :: a(*)
    integer(kind=8) :: i1x
!    CALCUL DE ALPHA*CONJG(Y)'
!-----------------------------------------------------------------------
! IN  : M    : NOMBRE DE LIGNE DE A.
!     : N    : NOMBRE DE COLONNE DE A.
!     : ALPHA: COMPLEXE.
!     : X    : VECTEUR COMPLEXE DE LONGUEUR (M-1)*IABS(INCX)+1.
!     : INCX : DEPLACEMENT ENTRE LES ELEMENTS DE X.
!     : Y    : VECTEUR COMPLEXE DE LONGUEUR (N-1)*IABS(INCY)+1.
!     : INCY : DEPLACEMENT ENTRE LES ELEMENTS DE Y.
! I/O : A    : MATRICE COMPLEXE DE DIMENSION M*N.
! IN  : LDA  : DIMENSION DE A.
!-----------------------------------------------------------------------
    integer(kind=8) :: iy, j
    blas_int :: b_incx, b_incy, b_n
    if (m .eq. 0 .or. n .eq. 0 .or. alpha .eq. (0.0d0, 0.0d0)) goto 9000
!
    iy = 1
    if (incy .lt. 0) iy = (-n+1)*incy+1
!
    i1x = 1
    do j = 1, n
        b_n = to_blas_int(m)
        b_incx = to_blas_int(incx)
        b_incy = to_blas_int(1)
        call zaxpy(b_n, alpha*dconjg(y(iy)), x, b_incx, a(i1x), &
                   b_incy)
        iy = iy+incy
        i1x = i1x+lda
    end do
!
9000 continue
    goto 999
999 continue
end subroutine
