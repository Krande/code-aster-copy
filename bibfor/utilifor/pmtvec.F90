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
subroutine pmtvec(cumul, n, a, x, y)
    implicit none
    character(len=*) :: cumul
    integer(kind=8) :: n
    real(kind=8) :: a(n, n), x(n), y(n)
!       ----------------------------------------------------------------
!       PRODUIT TRANSPOSEE D'UNE MATRICE CARREE PLEINE PAR UN VECTEUR
!          Y(N) = 0.  + A_t(N,N)*X(N)
!       OU Y(N) = Y(N)+ A_t(N,N)*X(N)
!       ----------------------------------------------------------------
! IN    N     : I  :   DIMENSION DE LA MATRICE ET DES VECTEURS X ET Y
! IN    A(N,N): R  :   MATRICE REELLE
! IN    X(N)  : R  :   VECTEUR REEL
! IN    CUMUL : K* :   ON CUMULE OU NON DANS LE VECTEUR RESULTAT Y
!       CUMUL = 'ZERO' ON MET Y A ZERO AVANT DE COMMENCER
!       CUMUL = 'CUMU' ON ACCUMULE DANS Y
! OUT   Y(N)  : R  :   VECTEUR REEL
!       ----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
!-----------------------------------------------------------------------
    if (cumul .eq. 'ZERO') then
        do i = 1, n
            y(i) = 0.d0
        end do
    end if
    do j = 1, n
        do i = 1, n
            y(i) = y(i)+a(j, i)*x(j)
        end do
    end do
end subroutine
