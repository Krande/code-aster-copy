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
subroutine mlmatc(ni, nk, nj, a, b, &
                  c)
    implicit none
    integer(kind=8) :: ni, nk, nj
    complex(kind=8) :: a(ni, *), b(nk, *), c(ni, *)
!     CALCUL COMPLEXE MATRICIEL C = A * B
!     MATRICE ORDONNEES PAR COLONNES DESCENDANTES
! ----------------------------------------------------------------------
! IN  : A  : MATRICE A(NI,NK)
!     : B  : MATRICE B(NK,NJ)
!     : C  : MATRICE C(NI,NJ)
!     : NI , NJ ,NK : DIMENSIONS DES MATRICES
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    complex(kind=8) :: xcres
!
    do i = 1, ni
        do j = 1, nj
            xcres = dcmplx(0.d0, 0.d0)
            do k = 1, nk
                xcres = xcres+a(i, k)*b(k, j)
            end do
            c(i, j) = xcres
        end do
    end do
!
end subroutine
