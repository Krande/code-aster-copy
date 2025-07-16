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
subroutine mltfld(n, front, adper, t1, ad, &
                  eps, ier)
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
#include "asterfort/sspmvb.h"
#include "blas/dgemv.h"
    integer(kind=8) :: n, adper(*), ad(*), ier
    real(kind=8) :: front(*), t1(*), eps, alpha, beta
    integer(kind=8) :: i, k
    integer(kind=8) :: seuin, seuik
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    parameter(seuin=1500, seuik=300)
    integer(kind=8) :: nn, kk, lda, incx, incy
    character(len=1) :: tra
    blas_int :: b_incx, b_incy, b_lda, b_m, b_n
!
!     MLTFLD TRAITE UNIQUEMENT LES BLOCS DIAGONAUX
!
    lda = n
    tra = 'N'
    alpha = -1.d0
    beta = 1.d0
    incx = 1
    incy = 1
    do k = 1, n
        do i = 1, k-1
            ad(i) = adper(i)+k-i
            t1(i) = front(ad(i))*front(adper(i))
        end do
        if (k .gt. 1) then
!
            nn = n-k+1
            kk = k-1
            if (nn .lt. seuin .or. kk .lt. seuik) then
                call sspmvb(n-k+1, k-1, front, ad, t1, &
                            front(adper(k)))
            else
                b_lda = to_blas_int(lda)
                b_m = to_blas_int(nn)
                b_n = to_blas_int(kk)
                b_incx = to_blas_int(incx)
                b_incy = to_blas_int(incy)
                call dgemv(tra, b_m, b_n, alpha, front(k), &
                           b_lda, t1, b_incx, beta, front(adper(k)), &
                           b_incy)
            end if
        end if
!         DIVISION PAR LE TERME DIAGONAL
        if (abs(front(adper(k))) .le. eps) then
            ier = k
            goto 40
        end if
        do i = 1, n-k
            front(adper(k)+i) = front(adper(k)+i)/front(adper(k))
        end do
    end do
40  continue
end subroutine
