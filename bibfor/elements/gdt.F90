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
subroutine gdt(teta, amat)
!
! FONCTION: LE VECTEUR-ROTATION ENTRE L'INSTANT INSTAM ET L'ITERATION
!           K-1 DE L'INSTANT INSTAP EST TETA. ON DONNE UN INCREMENT DE
!           VECTEUR-ROTATION DTETA A PARTIR DE L'ITERATION PRECEDENTE.
!           GDT CALCULE LA MATRICE AMAT QUI, MULTIPLIEE PAR DTETA, DONNE
!           LA PARTIE PRINCIPALE DE L'ACCROISSEMENT CORRESPONDANT DE
!           TETA.
!
!     IN  : TETA      : VECTEUR-ROTATION
!
!     OUT : AMAT      : MATRICE DE LINEARISATION
! ------------------------------------------------------------------
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/antisy.h"
#include "blas/ddot.h"
    real(kind=8) :: teta(3), eu(3), amat(3, 3), amat1(3, 3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
    real(kind=8) :: anor, anors2, coef, deux, epsil1, epsil2
    real(kind=8) :: prosca, un, zero
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    zero = 0.d0
    epsil1 = 1.d-4
    epsil2 = r8prem()**2
    un = 1.d0
    deux = 2.d0
!
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    prosca = ddot(b_n, teta, b_incx, teta, b_incy)
    anor = sqrt(prosca)
    anors2 = anor/deux
    if (anors2 .gt. epsil1) then
        coef = anors2/tan(anors2)
    else
        coef = un+anors2**2/3.d0
        coef = un/coef
    end if
!
    do i = 1, 3
        if (anors2 .le. epsil2) then
            eu(i) = zero
        else
            eu(i) = teta(i)/anor
        end if
    end do
    call antisy(teta, un, amat1)
    do j = 1, 3
        do i = 1, 3
            amat(i, j) = eu(i)*eu(j)*(un-coef)-amat1(i, j)/deux
        end do
        amat(j, j) = amat(j, j)+coef
    end do
end subroutine
