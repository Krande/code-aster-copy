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
subroutine fext(t, neq, nvect, liad, lifo, &
                f)
    implicit none
#include "jeveux.h"
#include "asterfort/fointe.h"
#include "asterfort/r8inir.h"
#include "blas/daxpy.h"
    integer(kind=8) :: neq, nvect, liad(*)
    real(kind=8) :: t, f(*)
    character(len=24) :: lifo(*)
!
!  CALCUL DU VECTEUR FEXT: FEXT = SOMME  GI(T)*FI(X)
!
!  INPUT:
!        T        : INSTANT DE CALCUL
!        NEQ      : NOMBRE D'EQUATIONS (D.D.L. ACTIFS)
!        NVECT    : NOMBRE DE VECTEURS CHARGEMENT
!        LIAD     : LISTE DES ADRESSES DES VECTEURS CHARGEMENT (NVECT)
!        LIFO     : LISTE DES NOMS DES FONCTIONS EVOLUTION (NVECT)
!
!  OUTPUT:
!        F        : VECTEUR FORCE EXTERIEURE (NEQ)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ier
    real(kind=8) :: zero, alpha
    character(len=8) :: nompar
    blas_int :: b_incx, b_incy, b_n
!     ------------------------------------------------------------------
!
    nompar = 'INST'
    zero = 0.d0
    call r8inir(neq, zero, f, 1)
    do i = 1, nvect
        call fointe('F ', lifo(i), 1, [nompar], [t], &
                    alpha, ier)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, alpha, zr(liad(i)), b_incx, f, &
                   b_incy)
    end do
end subroutine
