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
subroutine vroqua(theta, quater)
!
!
    implicit none
#include "blas/ddot.h"
    real(kind=8) :: quater(4), theta(3)
!
! ----------------------------------------------------------------------
!
! ROUTINE UTILITAIRE (QUATERNION)
!
! CALCULE LE QUATERNION CORRESPONDANT AU VECTEUR-ROTATION THETA
!
! ----------------------------------------------------------------------
!
!
! OUT QUATER : QUATERNION
! IN  THETA  : VECTEUR ROTATION
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: epsil, demi, un
    real(kind=8) :: angle, coef, prosca
    integer(kind=8) :: i
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
    epsil = 1.d-4
    demi = 5.d-1
    un = 1.d0
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    prosca = ddot(b_n, theta, b_incx, theta, b_incy)
    angle = demi*sqrt(prosca)
    quater(4) = cos(angle)
    if (angle .gt. epsil) then
        coef = demi*sin(angle)/angle
    else
        coef = un-angle**2/6.d0
        coef = demi*coef
    end if
    do i = 1, 3
        quater(i) = coef*theta(i)
    end do
!
end subroutine
