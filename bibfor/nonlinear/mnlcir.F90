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
subroutine mnlcir(xdep, ydep, omega, alpha, eta, &
                  h, hf, nt, xsort)
    implicit none
!
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnlfft.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
!
    integer(kind=8) :: h, hf, nt
    real(kind=8) :: omega, alpha, eta
    character(len=14) :: xdep, ydep, xsort
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    real(kind=8) :: depi, xt, yt, rk
    integer(kind=8) :: k, j, ix, iy, isor
    real(kind=8), pointer :: fn(:) => null()
    real(kind=8), pointer :: fx(:) => null()
    real(kind=8), pointer :: fy(:) => null()
    real(kind=8), pointer :: r(:) => null()
    real(kind=8), pointer :: t(:) => null()
    blas_int :: b_incx, b_n
!
    call jemarq()
!
    AS_ALLOCATE(vr=t, size=nt)
    AS_ALLOCATE(vr=fx, size=nt)
    AS_ALLOCATE(vr=fy, size=nt)
    AS_ALLOCATE(vr=fn, size=nt)
    AS_ALLOCATE(vr=r, size=nt)
!
    call jeveuo(xdep, 'E', ix)
    call jeveuo(ydep, 'E', iy)
    call jeveuo(xsort, 'E', isor)
!
    depi = r8depi()
    t(1) = 0.d0
    do k = 2, nt
        t(k) = t(k-1)+(depi/omega)/nt
    end do
!
    do k = 1, nt
        xt = zr(ix)
        yt = zr(iy)
        do j = 1, h
            xt = xt+zr(ix+j)*dcos(dble(j)*omega*t(k))
            yt = yt+zr(iy+j)*dcos(dble(j)*omega*t(k))
            xt = xt+zr(ix+h+j)*dsin(dble(j)*omega*t(k))
            yt = yt+zr(iy+h+j)*dsin(dble(j)*omega*t(k))
        end do
        rk = xt**2+yt**2
        r(k) = sqrt(rk)
        fn(k) = ((r(k)-1.d0)+sqrt((r(k)-1.d0)**2+4.d0*eta/alpha))/(2.d0/alpha)
        fx(k) = fn(k)*xt/r(k)
        fy(k) = fn(k)*yt/r(k)
    end do
!
    b_n = to_blas_int(4*(2*hf+1))
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(isor), b_incx)
    call mnlfft(1, zr(isor), fx, hf, nt, &
                1)
    call mnlfft(1, zr(isor+(2*hf+1)), fy, hf, nt, &
                1)
    call mnlfft(1, zr(isor+2*(2*hf+1)), r, hf, nt, &
                1)
    call mnlfft(1, zr(isor+3*(2*hf+1)), fn, hf, nt, &
                1)
!
    AS_DEALLOCATE(vr=t)
    AS_DEALLOCATE(vr=fx)
    AS_DEALLOCATE(vr=fy)
    AS_DEALLOCATE(vr=fn)
    AS_DEALLOCATE(vr=r)
!
    call jedema()
!
end subroutine
