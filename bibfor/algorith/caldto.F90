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
subroutine caldto(s6, fkooh, msns, dtods)
    implicit none
!
! person_in_charge: jean-michel.proix at edf.fr
!     ----------------------------------------------------------------
!
!     MONOCRISTAL : calcul de la derivee de Tau en GDEF : dTau/dS_ab
!     IN  S6    :  CONTRAINTES NOTATION VOIGT
!     IN  FKOOH :  INVERSE TENSEUR HOOKE
!     IN  MSNS  :  MS * NS
!     OUT DTODS :  dTau/dS
!
#include "asterfort/r8inir.h"
#include "asterfort/tnsvec.h"
#include "blas/dcopy.h"
    integer(kind=8) :: i, j, k, l, ind(3, 3), a, b, m, n
    real(kind=8) :: s6(6), fkooh(6, 6), msns(3, 3), dtods(3, 3)
    real(kind=8) :: s(3, 3), l4(3, 3, 3, 3)
    real(kind=8) :: mus(3, 3), l4s(3, 3)
    blas_int :: b_incx, b_incy, b_n
    data ind/1, 4, 5, 4, 2, 6, 5, 6, 3/
!     ----------------------------------------------------------------
!     CONSTUCTION DU TENSEUR INVERSE DE HOOKE D'ORDRE 4
!
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    l4(i, j, k, l) = fkooh(ind(i, j), ind(k, l))
                end do
            end do
        end do
    end do
!
    call tnsvec(6, 3, s, s6, 1.d0)
    b_n = to_blas_int(9)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, msns, b_incx, dtods, b_incy)
!
    call r8inir(9, 0.d0, mus, 1)
!
!     CALCUL DU TERME  2(Lambda**-1)*mu*S
    do m = 1, 3
        do n = 1, 3
            do k = 1, 3
                mus(m, n) = mus(m, n)+msns(m, k)*s(k, n)
            end do
        end do
    end do
!
    do a = 1, 3
        do b = 1, 3
            do m = 1, 3
                do n = 1, 3
                    dtods(a, b) = dtods(a, b)+2.d0*l4(a, b, m, n)*mus(m, n)
                end do
            end do
        end do
    end do
!
!     CALCUL DU TERME  2(LAMBDA**-1*S)*MU
    call r8inir(9, 0.d0, l4s, 1)
    do a = 1, 3
        do k = 1, 3
            do m = 1, 3
                do n = 1, 3
                    l4s(a, k) = l4s(a, k)+l4(a, k, m, n)*s(m, n)
                end do
            end do
        end do
    end do
!
    do a = 1, 3
        do b = 1, 3
            do k = 1, 3
                dtods(a, b) = dtods(a, b)+2.d0*l4s(a, k)*mus(k, b)
            end do
        end do
    end do
!
end subroutine
