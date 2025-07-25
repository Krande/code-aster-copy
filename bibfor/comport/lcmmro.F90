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
subroutine lcmmro(BEHinteg, omp, nvi, vind, vinf)
!
    use Behaviour_type
!
    implicit none
!
#include "asterc/r8miem.h"
#include "asterfort/r8inir.h"
#include "blas/dcopy.h"
!
!     Stockage variables internes rotation reseau
!     ----------------------------------------------------------------
    type(Behaviour_Integ), intent(in) :: BEHinteg
    integer(kind=8) :: i, j, nvi, k
    real(kind=8) :: omp(3), dtheta, iden(3, 3), nax(3, 3), q(3, 3)
    real(kind=8) :: omegap(3, 3), omegae(3, 3), omega(3, 3), dq(3, 3)
    real(kind=8) :: vind(nvi), vinf(nvi), l(3, 3), qm(3, 3)
    blas_int :: b_incx, b_incy, b_n
    data iden/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!
!     LA MATRICE DE ROTATION QM EST STOCKEE DANS VIND (N-19 A N-9)
    do i = 1, 3
        do j = 1, 3
            qm(i, j) = vind(nvi-19+3*(i-1)+j)+iden(i, j)
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            l(i, j) = BEHinteg%behavESVA%behavESVAGeom%gradVelo(3*(i-1)+j)
        end do
    end do
    do i = 1, 3
        do j = 1, 3
            omega(i, j) = 0.5d0*(l(i, j)-l(j, i))
        end do
    end do
!     LE VECTEUR  ROTATION PLASTIQUE EST STOCKE DANS VINF (N-9 A N-7)
    call r8inir(9, 0.d0, omegap, 1)
    omegap(2, 3) = -omp(1)
    omegap(3, 2) = +omp(1)
    omegap(1, 3) = +omp(2)
    omegap(3, 1) = -omp(2)
    omegap(1, 2) = -omp(3)
    omegap(2, 1) = +omp(3)
    do i = 1, 3
        do j = 1, 3
            omegae(i, j) = omega(i, j)-omegap(i, j)
        end do
    end do
!     ANGLE = NORME DU VECTEUR AXIAL
    dtheta = sqrt(omegae(1, 2)**2+omegae(1, 3)**2+omegae(2, 3)**2)
!
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, iden, b_incx, dq, b_incy)
    if (dtheta .gt. r8miem()) then
        do i = 1, 3
            do j = 1, 3
                nax(i, j) = omegae(i, j)/dtheta
            end do
        end do
        do i = 1, 3
            do j = 1, 3
                dq(i, j) = dq(i, j)+sin(dtheta)*nax(i, j)
            end do
        end do
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    dq(i, j) = dq(i, j)+(1.d0-cos(dtheta))*nax(i, k)*nax(k, j)
                end do
            end do
        end do
    end if
    call r8inir(9, 0.d0, q, 1)
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                q(i, j) = q(i, j)+dq(i, k)*qm(k, j)
            end do
        end do
    end do
!
! LA MATRICE DE ROTATION EST STOCKEE DANS VINF (N-18 A N-10)
    do i = 1, 3
        do j = 1, 3
            vinf(nvi-19+3*(i-1)+j) = (q(i, j)-iden(i, j))
        end do
    end do
!
! LE VECTEUR D-ROTATION PLASTIQUE EST STOCKE DANS VINF (N-9 A N-7)
!
    vinf(nvi-9) = omp(1)
    vinf(nvi-8) = omp(2)
    vinf(nvi-7) = omp(3)
!
! LE VECTEUR D-ROTATION ELASTIQUE EST STOCKE DANS VINF (N-6 A N-4)
!
    vinf(nvi-6) = omegae(3, 2)
    vinf(nvi-5) = omegae(1, 3)
    vinf(nvi-4) = omegae(2, 1)
!
    vinf(nvi-3) = dtheta+vind(nvi-3)
!
end subroutine
