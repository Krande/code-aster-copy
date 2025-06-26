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
subroutine caldfe(df, nr, nvi, vind, dfpds, &
                  fe, dfpdbs, msdgdt, drdy)
!
    implicit none
! person_in_charge: jean-michel.proix at edf.fr
!     ----------------------------------------------------------------
!
!     MONOCRISTAL : calcul des derivees de Fe en GDEF
!     IN  DF     :  GRADIENT DF
!         NR     :  DIMENSION DECLAREE DRDY
!         VIND   :  VARIABLES INTERNES A L'INSTANT PRECEDENT
!         DFPDS  :  DERIVEE DE FP PAR RAPPORT A S
!         FE     :  GRADIENT ELASTIQUE FE
!         DFPDBS :  DERIVEE DE FP PAR RAPPORT A BETA_S
!         MSDGDT :  SOMME DES MUS(I)*MUS(J)*DGAMMAS/DTAUS
!       OUT DRDY :  BLOC ((1-6),(7-NS)) JACOBIEN DU SYSTEME NON LINEAIRE
!
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    integer(kind=8) :: nr, ndt, ndi, ns, i, j, k, l, m, ind(3, 3), nvi
    real(kind=8) :: fe(3, 3), df(3, 3), dfpds(3, 3, 3, 3), msdgdt(6, 6)
    real(kind=8) :: dfefdt(3, 3, 3, 3)
    real(kind=8) :: vind(*), dfeds(3, 3, 3, 3), dfefds(3, 3, 3, 3), dffe(3, 3)
    real(kind=8) :: fem(3, 3)
    real(kind=8) :: id(3, 3), drdy(nr, nr)
    real(kind=8) :: dfpdbs(3, 3, 30), dfedbs(3, 3, 30), dfefdb(3, 3, 30)
    blas_int :: b_incx, b_incy, b_n
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
    data id/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
    data ind/1, 4, 5, 4, 2, 6, 5, 6, 3/
!     ----------------------------------------------------------------
!
    ns = nr-6
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
    dffe = matmul(df, fem)
!
!     on calcule dFe/dS
    call r8inir(81, 0.d0, dfeds, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        dfeds(i, j, k, l) = dfeds(i, j, k, l)+dffe(i, m)*dfpds(m, j, k, l)
                    end do
                end do
            end do
        end do
    end do
!
    call r8inir(81, 0.d0, dfefds, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        dfefds(i, j, k, l) = dfefds(i, j, k, l)+dfeds(m, i, k, l)*fe(m, j)
                    end do
                end do
            end do
        end do
    end do
!
    call r8inir(81, 0.d0, dfefdt, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    do m = 1, 3
                        dfefdt(i, j, k, l) = dfefdt(i, j, k, l)+dfeds(m, j, k, l)*fe(m, i)
                    end do
                end do
            end do
        end do
    end do
!
    b_n = to_blas_int(81)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 1.d0, dfefds, b_incx, dfefdt, &
               b_incy)
    b_n = to_blas_int(81)
    b_incx = to_blas_int(1)
    call dscal(b_n, -0.5d0, dfefdt, b_incx)
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    msdgdt(ind(i, j), ind(k, l)) = dfefdt(i, j, k, l)
                end do
            end do
        end do
    end do
!
!
!     on calcule dFe/dbetas
    call r8inir(3*3*ns, 0.d0, dfedbs, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 1, ns
                do m = 1, 3
                    dfedbs(i, j, k) = dfedbs(i, j, k)+dffe(i, m)*dfpdbs(m, j, k)
                end do
            end do
        end do
    end do
!
    call r8inir(3*3*ns, 0.d0, dfefdb, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 1, ns
                do m = 1, 3
                    dfefdb(i, j, k) = dfefdb(i, j, k)+dfedbs(m, i, k)*fe(m, j)+dfedbs(m, j, k)*fe&
                                      &(m, i)
                end do
            end do
        end do
    end do
!
!
    b_n = to_blas_int(3*3*ns)
    b_incx = to_blas_int(1)
    call dscal(b_n, -0.5d0, dfefdb, b_incx)
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, ns
                drdy(ind(i, j), 6+k) = dfefdb(i, j, k)
            end do
        end do
    end do
!
end subroutine caldfe
