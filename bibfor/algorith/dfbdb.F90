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
subroutine dfbdb(dim, b, e, deuxmu, lambda, &
                 ecrob, dsidep)
!
!
    implicit none
#include "asterfort/dfpdf.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: dim
    real(kind=8) :: b(6), e(6), deuxmu, lambda, dsidep(6, 6), ecrob
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT ENDO_ORTH_BETON
!     CALCUL DE LA DERIVEE DE LA FORCE THERMODYNAMIQUE (ENDO TRACTION)
!     PAR RAPPORT A L ENDOMMAGEMENT DE TRACTION:DFB/DB
!
!     FB=-LAMBDA.TR(EB).H(TR(EB))-MU/2*((BE+EB)_+*E + E*(BE+EB)_+)
!        +ECROB*(I-B)
!
!     IN  DIM      : DIMENSION 3(3D) OU 2(2D)
!     IN  E        : DEFORMATION
!     IN  B        : TENSEUR D ENDOMMAGEMENT DE TRACTION
!     IN  LAMBDA   : /
!     IN  DEUXMU   : / COEFFICIENTS DE LAME
!     IN  ECROB    : PARAMETRE D ECROUISSAGE DE TRACTION
!     OUT DSIDEP   : DFB/DB
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, l, t(3, 3), ndim
    real(kind=8) :: rtemp2
    real(kind=8) :: rtemp3, rtemp4, rac2, kron(3, 3), dbedb(6, 6), c(6)
    real(kind=8) :: mtemp(6, 6)
    real(kind=8) :: a(6, 6), treb
!
!
!
!
    if (dim .eq. 3) then
        t(1, 1) = 1
        t(1, 2) = 4
        t(1, 3) = 5
        t(2, 1) = 4
        t(2, 2) = 2
        t(2, 3) = 6
        t(3, 1) = 5
        t(3, 2) = 6
        t(3, 3) = 3
        ndim = 6
    else if (dim .eq. 2) then
        t(1, 1) = 1
        t(1, 2) = 3
        t(2, 2) = 2
        t(2, 1) = 3
        ndim = 3
    end if
!
    kron(1, 1) = 1.d0
    kron(1, 2) = 0.d0
    kron(1, 3) = 0.d0
    kron(2, 1) = 0.d0
    kron(2, 2) = 1.d0
    kron(2, 3) = 0.d0
    kron(3, 1) = 0.d0
    kron(3, 2) = 0.d0
    kron(3, 3) = 1.d0
!
    rac2 = sqrt(2.d0)
    call r8inir(36, 0.d0, dbedb, 1)
    call r8inir(6, 0.d0, c, 1)
    call r8inir(36, 0.d0, dsidep, 1)
    call r8inir(36, 0.d0, a, 1)
!
!
    do i = 1, dim
        do j = i, dim
            do k = 1, dim
                do l = 1, dim
                    if (i .eq. j) then
                        rtemp3 = 1.d0
                    else
                        rtemp3 = rac2
                    end if
                    if (k .eq. l) then
                        rtemp4 = 1.d0
                    else
                        rtemp4 = 1.d0/rac2
                    end if
                    dbedb(t(i, j), t(k, l)) = dbedb(t(i, j), t(k, l))+(kron(k, &
                                                  i)*e(t(l, j))+kron(j, l)*e(t(i, k)))*rtemp3*rtemp4
                end do
                c(t(i, j)) = c(t(i, j))+b(t(i, k))*e(t(k, j))+e(t(i, k))*b( &
                             t(k, j))
            end do
        end do
    end do
!
    call dfpdf(6, c, mtemp)
!
    do i = 1, ndim
        do j = 1, ndim
            do k = 1, ndim
                a(i, j) = a(i, j)+mtemp(i, k)*dbedb(k, j)
            end do
        end do
    end do
!
!
    do i = 1, dim
        do j = i, dim
            do k = 1, dim
                do l = 1, ndim
                    if (i .eq. j) then
                        rtemp2 = 1.d0
                    else
                        rtemp2 = rac2
                    end if
                    if (k .eq. j) then
                        rtemp3 = 1.d0
                    else
                        rtemp3 = 1.d0/rac2
                    end if
                    if (k .eq. i) then
                        rtemp4 = 1.d0
                    else
                        rtemp4 = 1.d0/rac2
                    end if
                    dsidep(t(i, j), l) = dsidep(t(i, j), l)+(a(t(i, k), l)*e( &
                                                 t(k, j))*rtemp4+e(t(i, k))*a(t(k, j), l)*rtemp3)* &
                                         rtemp2
                end do
            end do
        end do
    end do
!
!
    do i = dim+1, ndim
        e(i) = rac2*e(i)
    end do
!
    if (dim .eq. 3) then
        treb = (c(1)+c(2)+c(3))
    else if (dim .eq. 2) then
        treb = (c(1)+c(2))
    end if
!
!
    if (treb .lt. 0.d0) then
        treb = 0.d0
    else
        treb = 1.d0
    end if
    do i = 1, ndim
        do j = 1, ndim
            dsidep(i, j) = -deuxmu/4.d0*dsidep(i, j)-lambda*treb*e(i)*e(j)
        end do
    end do
!
    do i = dim+1, ndim
        e(i) = e(i)/rac2
    end do
!
!
!CC ON RAJOUTE LE TERME QUI VIENT DE L ECROUISSAGE
!
    do i = 1, ndim
        dsidep(i, i) = dsidep(i, i)-ecrob
    end do
!
!
!
end subroutine
