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
subroutine zmvpy(uplo, n, alpha, a, lda, &
                 x, incx, beta, y, incy)
    implicit none
#include "asterfort/vecinc.h"
#include "asterfort/zmult.h"
#include "blas/zaxpy.h"
    integer(kind=8) :: n, lda, incx, incy
    complex(kind=8) :: alpha, beta, a(lda, *), x(*), y(*)
    character(len=*) :: uplo
!  PRODUIT D'UNE MATRICE HERMITIENNE PAR UN VECTEUR SUIVANT OPTION
!                         'U' OU 'L'.
!-----------------------------------------------------------------------
! IN  : UPLO : CARACTERE SPECIFIANT LE MODE DE STOCKAGE.
!                 UPLO              STRUCTURE
!              'U'             MATRICE TRIANGULAIRE SUPERIEURE.
!              'L'             MATRICE TRIANGULAIRE INFERIEURE
!     : N    : DIMENSION DE LA MATRICE.
!     : ALPHA: COMPLEXE.
!     : A    : MATRICE COMPLXE DE DIMENSION M*N.
!     : LDA  : DIMENSION DE A.
!     : X    : DVECTEUR COMPLXE DE LONGUEUR (N-1)*IABS(INCX)+1.
!     : INCX : DEPLACEMENT ENTRE LES ELEMENTS DE X.
!     : BETA : COMPLEXE.LORSQUE BETA EGAL ZERO, Y N'EST PAS CALCULE.
!     : Y    : VECTEUR COMPLEXE DE LONGUEUR (N-1)*IABS(INCY)+1.
!     : INCY : DEPLACEMENT ENTRE LES ELEMENTS DE Y.
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ix, iy, j, ky
    complex(kind=8) :: temp
    real(kind=8) :: dble
    blas_int :: b_incx, b_incy, b_n
!
    if (n .eq. 0 .or. (alpha .eq. (0.0d0, 0.0d0) .and. beta .eq. (1.0d0, 0.0d0))) goto 999
!
    ix = 1
    iy = 1
    if (incx .lt. 0) ix = (-n+1)*incx+1
    if (incy .lt. 0) iy = (-n+1)*incy+1
!
    if (beta .eq. (1.0d0, 0.0d0)) then
    else if (incy .eq. 0) then
        if (beta .eq. (0.0d0, 0.0d0)) then
            y(1) = (0.0d0, 0.0d0)
        else
            y(1) = beta**n*y(1)
        end if
    else if (beta .eq. (0.0d0, 0.0d0)) then
        call vecinc(n, (0.0d0, 0.0d0), y, inc=abs(incy))
    else
        call zmult(n, beta, y, abs(incy))
    end if
!
    if (alpha .eq. (0.0d0, 0.0d0)) goto 999
!
    if (uplo(1:1) .eq. 'U' .or. uplo(1:1) .eq. 'u') then
        do j = 1, n
            temp = alpha*x(ix)
            ky = iy+(j-2)*min(incy, 0)
            b_n = to_blas_int(j-1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(incy)
            call zaxpy(b_n, temp, a(1, j), b_incx, y(ky), &
                       b_incy)
            ky = iy+(j-1)*incy
            y(ky) = y(ky)+temp*dble(a(j, j))
            do i = j+1, n
                ky = ky+incy
                y(ky) = y(ky)+temp*dconjg(a(j, i))
            end do
            ix = ix+incx
        end do
    else
        do j = 1, n
            temp = alpha*x(ix)
            ky = iy
            do i = 1, j-1
                y(ky) = y(ky)+temp*dconjg(a(j, i))
                ky = ky+incy
            end do
            y(ky) = y(ky)+temp*dble(a(j, j))
            ky = ky+incy+(n-j-1)*min(incy, 0)
            b_n = to_blas_int(n-j)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(incy)
            call zaxpy(b_n, temp, a(j+1, j), b_incx, y(ky), &
                       b_incy)
            ix = ix+incx
        end do
    end if
!
999 continue
end subroutine
