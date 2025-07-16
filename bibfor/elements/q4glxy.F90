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
subroutine q4glxy(hlt2, depf, lambda)
    implicit none
    real(kind=8) :: hlt2(4, 6), depf(12), lambda(4)
!     'LAMBDA' DE L'ELEMENT DE PLAQUE Q4G
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    real(kind=8) :: tb(6, 12)
    real(kind=8) :: blb(4, 12)
!     ------------------------------------------------------------------
!
!       ---- CALCUL DE LA MATRICE TB -------------------------------
    do k = 1, 6
        do j = 1, 12
            tb(k, j) = 0.d0
        end do
    end do
    tb(3, 2) = 0.25d0
    tb(3, 5) = -0.25d0
    tb(3, 8) = 0.25d0
    tb(3, 11) = -0.25d0
    tb(6, 3) = 0.25d0
    tb(6, 6) = -0.25d0
    tb(6, 9) = 0.25d0
    tb(6, 12) = -0.25d0
!
!        -------------- BLB = HLT2.TB ---------------------------
    do i = 1, 4
        do j = 1, 12
            blb(i, j) = 0.d0
            do k = 1, 6
                blb(i, j) = blb(i, j)+hlt2(i, k)*tb(k, j)
            end do
        end do
    end do
!        -------- LAMBDA = BLB.DEPF -----------------------------
    do i = 1, 4
        lambda(i) = 0.d0
    end do
    do i = 1, 4
        do j = 1, 12
            lambda(i) = lambda(i)+blb(i, j)*depf(j)
        end do
    end do
!
end subroutine
