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
subroutine mtcrou(a, b, nmax, n, nbscmb, &
                  l, d)
    implicit none
    integer(kind=8) :: n, nbscmb, nmax
    real(kind=8) :: a(nmax, nmax), b(nmax, nbscmb), l(n, n), d(n)
!     ------------------------------------------------------------------
!     RESOLUTION PAR LA METHODE DE CROUT D'UN SYSTEME LINEAIRE
!     ------------------------------------------------------------------
! VAR A      : R8 : MATRICE CARREE PLEINE
! VAR B      : R8 : TABLEAU BI-INDICES DE REELS
!               EN ENTREE : LES SECONDS MEMBRES
!               EN SORTIE : LES SOLUTIONS
! IN  NMAX   : IS : DIM MAXI DE LA MATRICE
! IN  N      : IS : ORDRE DE LA MATRICE
! IN  NBSCMB : IS : NOMBRE DE SECOND MEMBRE
!     ------------------------------------------------------------------
    real(kind=8) :: zero, s
    integer(kind=8) :: i, is, j, k
!-----------------------------------------------------------------------
    zero = 0.d0
    do i = 1, n
        do j = 1, i-1
            s = zero
            do k = 1, j-1
                s = s+l(i, k)*d(k)*l(j, k)
            end do
            l(i, j) = (a(i, j)-s)/d(j)
        end do
        s = zero
        do k = 1, i-1
            s = s+l(i, k)*l(i, k)*d(k)
        end do
        d(i) = a(i, i)-s
    end do
!
!   BOUCLE SUR LES SECONDS MEMBRES
!
    do is = 1, nbscmb
!
!  DESCENTE
!
        do i = 1, n
            s = zero
            do k = 1, i-1
                s = s+l(i, k)*b(k, is)
            end do
            b(i, is) = b(i, is)-s
        end do
!
!  DIVISION PAR LA DIAGONALE
!
        do i = 1, n
            b(i, is) = b(i, is)/d(i)
        end do
!
!  REMONTEE
!
        do i = n, 1, -1
            s = zero
            do k = i+1, n
                s = s+l(k, i)*b(k, is)
            end do
            b(i, is) = b(i, is)-s
        end do
    end do
end subroutine
