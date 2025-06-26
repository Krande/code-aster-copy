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
subroutine forcdy(masse, amort, lamort, neq, c0, &
                  c1, c2, c3, c4, c5, &
                  d0, v0, a0, f1, f2, &
                  f)
!   BUT :      CALCUL DU VECTEUR FORCE DYNAMIQUE
!
!              F  = F  + M*(C0.D0+C1.V0+C2.A0)
!                      + C*(C3.D0+C4.V0+C5.A0)
! ----------------------------------------------------------------------
!   INPUT:
!   ---> MASSE   : POINTEUR DE LA MATRICE MASSE
!   ---> AMORT   : POINTEUR DE LA MATRICE AMORTISS
!   ---> LAMORT  : VARIABLE LOGIQUE
!                     .TRUE. SI IL Y A UNE MATRICE AMORTISSEMENT
!                     .FALSE. SINON
!   ---> C0,C1,C2,C3,C4,C5 : CONSTANTES DE CALCUL
!   ---> NEQ   : NOMBRE D'EQUATIONS (D.D.L. ACTIFS)
!   ---> D0    : VECTEUR DEPLACEMENT  INITIAL  (NEQ)
!   ---> V0    : VECTEUR VITESSE      INITIAL  (NEQ)
!   ---> A0    : VECTEUR ACCELERATION INITIAL  (NEQ)
!   ---> F1    : VECTEUR REEL DE TRAVAIL        (NEQ)
!   ---> F2    : VECTEUR REEL DE TRAVAIL        (NEQ)
!
!   VAR   :
!   <--> F     : VECTEUR FORCE EXTERIEURE ENTREE (NEQ)
!                VECTEUR FORCE DYNAMIQUE SORTIE (NEQ)
! ----------------------------------------------------------------------
!
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
!
#include "asterf_types.h"
#include "asterfort/mrmult.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
    integer(kind=8) :: masse, amort
    integer(kind=8) :: neq
!
    real(kind=8) :: d0(*), v0(*), a0(*), f1(*), f2(*), f(*)
    real(kind=8) :: c0, c1, c2, c3, c4, c5
!
    aster_logical :: lamort
!
! DECLARATION VARIABLES LOCALES
!
    real(kind=8) :: zero, un
    blas_int :: b_incx, b_incy, b_n
!
    zero = 0.d0
    un = 1.d0
!
! 1. TERME M*(C0.D0+C1.V0+C2.A0)
!
    call r8inir(neq, zero, f1, 1)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, c0, d0, b_incx, f1, &
               b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, c1, v0, b_incx, f1, &
               b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, c2, a0, b_incx, f1, &
               b_incy)
!
    call mrmult('ZERO', masse, f1, f2, 1, &
                .true._1)
!
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, un, f2, b_incx, f, &
               b_incy)
!
! 2. CUMUL EVENTUEL DE C*(C3.D0+C4.V0+C5.A0)
!
    if (lamort) then
!
        call r8inir(neq, zero, f1, 1)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, c3, d0, b_incx, f1, &
                   b_incy)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, c4, v0, b_incx, f1, &
                   b_incy)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, c5, a0, b_incx, f1, &
                   b_incy)
!
        call mrmult('ZERO', amort, f1, f2, 1, &
                    .true._1)
!
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, un, f2, b_incx, f, &
                   b_incy)
!
    end if
!
end subroutine
