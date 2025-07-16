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
subroutine quavro(quater, theta)
!
!
! aslint: disable=
    implicit none
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "blas/ddot.h"
    real(kind=8) :: quater(4), theta(3)
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (QUATERNION)
!
! CALCULE LE VECTEUR-ROTATION THETA CORRESPONDANT AU QUATERNION
! QUATER
!
! ----------------------------------------------------------------------
!
!
! IN  QUATER : QUATERNION
! OUT THETA  : VECTEUR ROTATION
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: reste, zero, epsil, deux, coef
    real(kind=8) :: pi
    real(kind=8) :: prosca, anorx
    integer(kind=8) :: i
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!C
    zero = 0.d0
    epsil = r8prem()**2
    deux = 2.d0
    pi = r8pi()
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    prosca = ddot(b_n, quater, b_incx, quater, b_incy)
    anorx = sqrt(prosca)
    if (anorx .gt. 1.d0) anorx = 1.d0
    if (anorx .lt. epsil) then
        do i = 1, 3
            theta(i) = zero
        end do
        goto 999
    end if
    reste = asin(anorx)
!
    if (quater(4) .lt. zero) reste = pi-reste
!
    coef = deux*reste/anorx
!
    do i = 1, 3
        theta(i) = coef*quater(i)
    end do
999 continue
end subroutine
