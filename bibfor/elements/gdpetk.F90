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
subroutine gdpetk(tetag, tetapg, petikm, petik)
!
! FONCTION: POUR UN ELEMENT DE POUTRE EN GRAND DEPLACEMENT, CALCULE LE
!           VECTEUR-COURBURE AUX POINTS DE GAUSS, EN AXES GENERAUX,
!           VECTEUR DENOMME 'PETIT K'.
!
!     IN  : TETAG     : VECTEUR-INCREMENT DE ROTATION A GAUCHE
!           TETAPG    : DERIVEE DE TETAG PAR RAPPORT A L'ABS. CURVILIGNE
!           PETIKM    : VECTEUR-COURBURE A L'ITERATION PRECEDENTE
!
!     OUT : PETIK     : VECTEUR-COURBURE ACTUEL
! ------------------------------------------------------------------
    implicit none
#include "asterfort/antisy.h"
#include "asterfort/axial.h"
#include "asterfort/marota.h"
#include "asterfort/promat.h"
#include "asterfort/provec.h"
#include "asterfort/transp.h"
#include "blas/ddot.h"
    real(kind=8) :: tetag(3), tetapg(3), petikm(3), petik(3), petik1(3)
    real(kind=8) :: petik2(3), v1(3), amat1(3, 3), amat2(3, 3), amat3(3, 3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i
    real(kind=8) :: coef1, coef2, coef3, demi, epsil, prosca
    real(kind=8) :: teta1, teta2, un
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    epsil = 1.d-8
    demi = 5.d-1
    un = 1.d0
!
!*** PETIK1: VECTEUR BETA (SIMO: 'A THREE-DIMENSIONAL FINITE-STRAIN ROD
!***                       MODEL-PART 2'.)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    teta2 = ddot(b_n, tetag, b_incx, tetag, b_incy)
    if (abs(teta2) .lt. epsil) goto 11
    teta1 = sqrt(teta2)
    call provec(tetag, tetapg, v1)
    coef1 = sin(teta1)/teta1
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    prosca = ddot(b_n, tetag, b_incx, tetapg, b_incy)
    coef2 = (un-coef1)*prosca/teta2
    coef3 = demi*(sin(demi*teta1)/(demi*teta1))**2
    do i = 1, 3
        petik1(i) = coef1*tetapg(i)+coef2*tetag(i)+coef3*v1(i)
    end do
    goto 20
!
!*** TETAG EST TRES PETIT ET BETA VAUT PRATIQUEMENT TETAPRIM
11  continue
    do i = 1, 3
        petik1(i) = tetapg(i)
    end do
!
20  continue
!
!
    call marota(tetag, amat1)
    call antisy(petikm, un, amat2)
    call promat(amat1, 3, 3, 3, amat2, &
                3, 3, 3, amat3)
    call transp(amat1, 3, 3, 3, amat2, &
                3)
    call promat(amat3, 3, 3, 3, amat2, &
                3, 3, 3, amat1)
    call axial(amat1, petik2)
!
    do i = 1, 3
        petik(i) = petik1(i)+petik2(i)
    end do
end subroutine
