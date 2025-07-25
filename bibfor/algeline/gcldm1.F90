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
subroutine gcldm1(m, in, ip, prec, x, &
                  y, perm, xtrav, ytrav)
    implicit none
!
    integer(kind=8), intent(in) :: m
    real(kind=8), intent(in) :: prec(*), x(*)
    real(kind=8), intent(out) :: y(*)
    integer(kind=4), intent(in) :: ip(*)
    integer(kind=8), intent(in) :: in(*)
    integer(kind=8), intent(in) :: perm(*)
    real(kind=8), intent(inout) :: xtrav(*)
    real(kind=8), intent(inout) :: ytrav(*)
!-----------------------------------------------------------------------
!  FONCTION  :  INVERSION D'UNE MATRICE DE PRECONDITIONNEMENT LDLT_INC
!                        -1                   T
!               Y = (MAT)  *X   OU MAT = L*D*L
!          LA MATRICE MAT EST STOCKEE SOUS FORME MORSE
!                                                      -1
!          ET A LA PLACE DE 1. DANS MAT(I,I) ON A (D(I))
!               RESOLUTION DU SYSTEME :
!                                    T -1
!                          Y = (L D L ) * X
!-----------------------------------------------------------------------
    integer(kind=8) :: i, kdeb, kfin, ki
    real(kind=8) :: som, fac
!-----------------------------------------------------------------------
!
!   -- on permute x pour qu'il ait la numerotation du preconditionneur :
!   --------------------------------------------------------------------
    do i = 1, m
        xtrav(perm(i)) = x(i)
    end do
!
!-----------------------------------------------------------------------
!     RESOLUTION DU PREMIER SYSTEME L.W = X
!-------------------------------------------
    ytrav(1) = xtrav(1)
!
    do i = 2, m
        som = 0.d0
        kdeb = in(i-1)+1
        kfin = in(i)-1
        do ki = kdeb, kfin
            som = som+prec(ki)*ytrav(ip(ki))
        end do
        ytrav(i) = (xtrav(i)-som)
    end do
!
!-------------------------------------------
!     RESOLUTION DE D.Y = W
!-------------------------------------------
    do i = 1, m
        ytrav(i) = ytrav(i)*prec(in(i))
    end do
!
!-------------------------------------------
!     RESOLUTION DU SECOND SYSTEME LT.Y = W
!-------------------------------------------
    do i = m, 2, -1
        kdeb = in(i-1)+1
        kfin = in(i)-1
        fac = ytrav(i)
!
!        ---- PROCEDURE A LA MAIN
        do ki = kdeb, kfin
            ytrav(ip(ki)) = ytrav(ip(ki))-prec(ki)*fac
        end do
    end do
!
!   -- on permute ytrav pour qu'il ait la numerotation du syteme :
!   --------------------------------------------------------------
    do i = 1, m
        y(i) = ytrav(perm(i))
    end do
!
!
end subroutine
