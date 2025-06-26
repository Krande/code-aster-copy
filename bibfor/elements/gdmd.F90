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
subroutine gdmd(x0pg, pn, pm, d)
!
! FONCTION: POUR UN ELEMENT DE POUTRE EN GRAND DEPLACEMENT, CALCULE LA
!           MATRICE D AUX POINTS DE GAUSS. LA MATRICE D INTERVIENT DANS
!           LA RIGIDITE GEOMETRIQUE COMME LA MATRICE DE COMPORTEMENT C
!           INTERVIENT DANS LA RIGIDITE MATERIELLE.
!
!     IN  : VALEURS AU POINT DE GAUSS
!           X0PG      : DERIVEES DES COORDONNEES PAR RAP. A L'ABS. CURV.
!           PN        : RESULTANTE DES FORCES EN AX.GENE.
!           PM        : MOMENT RESULTANT EN AXES GENERAUX
!
!     OUT : D         : MATRICE 9*9
! ------------------------------------------------------------------
    implicit none
#include "asterfort/antisy.h"
#include "blas/ddot.h"
    real(kind=8) :: x0pg(3), pn(3), pm(3), d(9, 9), pntild(3, 3), pmtild(3, 3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
    real(kind=8) :: scal, un, zero
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    zero = 0.d0
    un = 1.d0
    do j = 1, 9
        do i = 1, 9
            d(i, j) = zero
        end do
    end do
    call antisy(pn, un, pntild)
    call antisy(pm, un, pmtild)
    do i = 1, 3
        do j = 1, 3
            d(i, 6+j) = -pntild(i, j)
            d(3+i, 6+j) = -pmtild(i, j)
            d(6+i, j) = pntild(i, j)
        end do
    end do
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    scal = ddot(b_n, pn, b_incx, x0pg, b_incy)
    do j = 1, 3
        do i = 1, 3
            d(6+i, 6+j) = pn(i)*x0pg(j)
        end do
        d(6+j, 6+j) = d(6+j, 6+j)-scal
    end do
end subroutine
