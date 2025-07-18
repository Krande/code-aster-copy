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
subroutine lglpma(n, a, b, c)
!
    implicit none
    integer(kind=8) :: n
    real(kind=8) :: a(6, 6), b(6, 6), c(6, 6)
! --- BUT : CALCUL DU PRODUIT MATRICIEL --------------------------------
! ======================================================================
! IN  : N      : DIMENSION REELLE DES MATRICES -------------------------
! --- : DIM    : DIMENSION EFFECTIVE DES MATRICES ----------------------
! --- : A      : MATRICE A ---------------------------------------------
! --- : B      : MATRICE B ---------------------------------------------
! OUT : C      : MATRICE C = A * B -------------------------------------
! ======================================================================
    integer(kind=8) :: i, j, k
! ======================================================================
    real(kind=8) :: v
    do i = 1, n
        do j = 1, n
            v = 0.0d0
            do k = 1, n
                v = v+a(i, k)*b(k, j)
            end do
            c(i, j) = v
        end do
    end do
! ======================================================================
end subroutine
