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
subroutine gcax(m, in, ip, ac, x, &
                y)
    implicit none
    integer(kind=4) :: ip(*)
    integer(kind=8) :: m, in(m)
    real(kind=8) :: ac(*), x(m), y(m)
    real(kind=8) :: dtemp
!     ------------------------------------------------------------------
!     MULTIPLICATION D'UNE MATRICE SYMETRIQUE COMPACTE PAR
!                UN VECTEUR :  Y = AC*X
!     ------------------------------------------------------------------
! IN . M             -->   NOMBRE DE COLONNES DE LA MATRICE
! IN . IN(I=1,M)     -->   POINTEUR DE FIN DE COLONNE DE LA MATRICE
! IN . IP(J)         -->   TABLEAU DES NUMEROS DE LIGNE
! IN . AC(J)         -->   TABLEAU DES COEFFICIENTS DE LA MATRICE
! IN . X(I=1,M)      -->   VECTEUR D'ENTREE
! OUT. Y(I=1,M)     <--    VECTEUR DE SORTIE
!     _____________ ____ ______________________________________________
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, kdeb, kfin, ki, klong
!-----------------------------------------------------------------------
    y(1) = ac(1)*x(1)
    do i = 2, m
        kdeb = in(i-1)+1
        kfin = in(i)-1
        klong = in(i)
        dtemp = 0.0d0
        do j = kdeb, klong
            dtemp = dtemp+x(ip(j))*ac(j)
        end do
        y(i) = dtemp
        dtemp = x(i)
        do ki = kdeb, kfin
            y(ip(ki)) = y(ip(ki))+ac(ki)*dtemp
        end do
    end do
!
end subroutine
