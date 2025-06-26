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
subroutine lcmmkg(zinv, nvi, vind, vinf, nmat, &
                  materf, mod, nr, dsde)
    implicit none
!
#include "asterfort/lcopli.h"
#include "asterfort/matinv.h"
#include "asterfort/r8inir.h"
#include "asterfort/tnsvec.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    real(kind=8) :: fep(3, 3), fem(3, 3), vinf(*), vind(*), dsde(6, 3, 3), det
    real(kind=8) :: dfedf(3, 3, 3, 3), fpp(3, 3), fppinv(3, 3), id(3, 3)
    real(kind=8) :: dr1df(3, 3, 3, 3), dr1df6(6, 3, 3), zinv(6, 6)
    real(kind=8) :: dsdf(6, 3, 3)
    real(kind=8) :: fet(3, 3), fetfe(3, 3), eel(6), hooke(6, 6), s6(6), s(3, 3)
    real(kind=8) :: dtaudf(3, 3, 3, 3), materf(*), sfet(3, 3)
    integer(kind=8) :: nr, ndt, ndi, i, j, k, l, m, n, nmat, ind(3, 3), nvi
    common/tdim/ndt, ndi
    character(len=8) :: mod
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
    data id/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!
    do i = 1, 3
        ind(i, i) = i
    end do
    ind(1, 2) = 4
    ind(2, 1) = 4
    ind(1, 3) = 5
    ind(3, 1) = 5
    ind(2, 3) = 6
    ind(3, 2) = 6
!
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vind(nvi-3-18+10), b_incx, fem, b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, id, b_incx, fem, &
               b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vinf(nvi-3-18+10), b_incx, fep, b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, id, b_incx, fep, &
               b_incy)
!
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vinf(nvi-3-18+1), b_incx, fpp, b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, id, b_incx, fpp, &
               b_incy)
    call matinv('S', 3, fpp, fppinv, det)
!
! CALCUL DE DFE/DF
!
    call r8inir(81, 0.d0, dfedf, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        dfedf(i, j, k, l) = dfedf(i, j, k, l)+id(i, k)*fem(l, m)*fppinv(m, j)
                    end do
                end do
            end do
        end do
    end do
!
! CALCUL DE DR1/DF
    call r8inir(81, 0.d0, dr1df, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        dr1df(i, j, k, l) = dr1df(i, j, k, l)+dfedf(m, i, k, l)*fep(m, j)
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        dr1df(i, j, k, l) = dr1df(i, j, k, l)+fep(m, i)*dfedf(m, j, k, l)
                    end do
                end do
            end do
        end do
    end do
    b_n = to_blas_int(81)
    b_incx = to_blas_int(1)
    call dscal(b_n, -0.5d0, dr1df, b_incx)
!
! CALCUL DE DS/DF EN UTILISANT LES SYMETRIES
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    dr1df6(ind(i, j), k, l) = dr1df(i, j, k, l)
                end do
            end do
        end do
    end do
!
    call r8inir(54, 0.d0, dsdf, 1)
    do i = 1, 6
        do j = 1, 3
            do k = 1, 3
                do l = 1, 6
                    dsdf(i, j, k) = dsdf(i, j, k)-zinv(i, l)*dr1df6(l, j, k)
                end do
            end do
        end do
    end do
!
! RECALCUL DU PK2 S
    if (materf(nmat) .eq. 0) then
        call lcopli('ISOTROPE', mod, materf(1), hooke)
    else if (materf(nmat) .eq. 1) then
        call lcopli('ORTHOTRO', mod, materf(1), hooke)
    end if
    fet = transpose(fep)
    fetfe = matmul(fet, fep)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -1.d0, id, b_incx, fetfe, &
               b_incy)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.5d0, fetfe, b_incx)
!
!      CONTRAINTES PK2
    call tnsvec(3, 3, fetfe, eel, 1.d0)
    s6(1:ndt) = matmul(hooke(1:ndt, 1:ndt), eel(1:ndt))
    call tnsvec(6, 3, s, s6, 1.d0)
!
! CALCUL DE DTAU/DF EN UTILISANT LES SYMETRIES
    call r8inir(81, 0.d0, dtaudf, 1)
    sfet = matmul(s, fet)
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        dtaudf(i, j, k, l) = dtaudf(i, j, k, l)+dfedf(i, m, k, l)*sfet(m, j)
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        dtaudf(i, j, k, l) = dtaudf(i, j, k, l)+sfet(m, i)*dfedf(j, m, k, l)
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        do n = 1, 3
                            dtaudf(i, j, k, l) = dtaudf(i, j, k, l)+fep(i, m)*dsdf(ind(m, n), k, &
                                                 &l)*fep(j, n)
                        end do
                    end do
                end do
            end do
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    dsde(ind(i, j), k, l) = dtaudf(i, j, k, l)
                end do
            end do
        end do
    end do
!
! LES RACINE(2) ATTENDUES PAR NMCOMP  !!!
    b_n = to_blas_int(9)
    b_incx = to_blas_int(6)
    call dscal(b_n, sqrt(2.d0), dsde(4, 1, 1), b_incx)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(6)
    call dscal(b_n, sqrt(2.d0), dsde(5, 1, 1), b_incx)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(6)
    call dscal(b_n, sqrt(2.d0), dsde(6, 1, 1), b_incx)
!
end subroutine
