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
subroutine newacc(neq, c1, c2, c3, d0, &
                  v0, a0, d1, a1)
    implicit none
!
!
!  INPUT:
!        NEQ   : NOMBRE D'EQUATIONS (D.D.L. ACTIFS)
!        C1,C2,C3 : CONSTANTES DE CALCUL
!        D0    : VECTEUR DEPLACEMENT  INITIAL   (NEQ)
!        V0    : VECTEUR VITESSE      INITIALE  (NEQ)
!        A0    : VECTEUR ACCELERATION INITIALE  (NEQ)
!        D1    : VECTEUR DEPLACEMENT SOLUTION   (NEQ)
!
!  OUTPUT:
!        A1    : VECTEUR ACCEL;ERATION SOLUTION (NEQ)
!
!  CALCUL DE L'ACCELERATION: ACCSOL = C1*( DEPSOL-D0) + C2*V0 + C3*A0
!
!
!----------------------------------------------------------------------
!   E.D.F DER   JACQUART G. 47-65-49-41      LE 19 JUILLET 1990
!**********************************************************************
!
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
    real(kind=8) :: d0(*), d1(*), v0(*), a0(*), a1(*)
    real(kind=8) :: c1, c2, c3
    real(kind=8) :: scal
!
!-----------------------------------------------------------------------
    integer(kind=8) :: neq
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    scal = -1.d0
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, d1, b_incx, a1, b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, scal, d0, b_incx, a1, &
               b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    call dscal(b_n, c1, a1, b_incx)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, c2, v0, b_incx, a1, &
               b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, c3, a0, b_incx, a1, &
               b_incy)
end subroutine
