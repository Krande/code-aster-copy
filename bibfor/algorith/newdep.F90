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
subroutine newdep(neq, c, dt, d0, v0, &
                  a0, d1, a1)
    implicit none
!
!
!   INPUT:
!        NEQ   : NOMBRE D'EQUATIONS (D.D.L. ACTIFS)
!        C     : CONSTANTE DE CALCUL
!        DT    : PAS DE TEMPS DE CALCUL
!        D0    : VECTEUR DEPLACEMENT  INITIAL  (NEQ)
!        V0    : VECTEUR VITESSE      INITIALE (NEQ)
!        A0    : VECTEUR ACCELERATION INITIALE (NEQ)
!        A1    : VECTEUR ACCELERATION SOLUTION (NEQ)
!
!   OUTPUT:
!        D1    : VECTEUR DEPLACEMENT  SOLUTION (NEQ)
!
!  CALCUL DU DEPLACEMENT (WILSON)
!  ============================= DEPSOL = D0 + DT*V0 + C*(ACCSOL+2.A0)
!
!
!----------------------------------------------------------------------
!   E.D.F DER   JACQUART G. 47-65-49-41      LE 19 JUILLET 1990
!**********************************************************************
!
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    real(kind=8) :: d0(*), d1(*), v0(*), a0(*), a1(*), c, dt
!-----------------------------------------------------------------------
    integer(kind=8) :: neq
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, d0, b_incx, d1, b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, dt, v0, b_incx, d1, &
               b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, c, a1, b_incx, d1, &
               b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, 2*c, a0, b_incx, d1, &
               b_incy)
end subroutine
