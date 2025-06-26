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
subroutine amppr(amat, nb1, nb2, bmat, n1, &
                 n2, i, j)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 12/03/91
!-----------------------------------------------------------------------
!  BUT:  AJOUTER UNE MATRICE PLEINE REELLE A UNE MATRICEPLEINE REELLE
!            A UNE MATRICE PLEINE REELLE
!-----------------------------------------------------------------------
!
! AMAT     /M/: MATRICE RECEPTRICE
! NB1      /I/: NB DE LIGNES DE LA MATRICE RECEPTRICE
! NB2      /I/: NB DE COLONNES DE LA MATRICE RECEPTRICE
! BMAT     /M/: MATRICE PLEINE A AJOUTER
! N1       /I/: NB DE LIGNE DE LA MATRICE A AJOUTER
! N2       /I/: NB DE COLONNE DE LA MATRICE A AJOUTER
! I        /I/: INDICE DU PREMIER TERME DANS RECEPTRICE
! J        /I/: INDICE DE COLONNE TERME  DANS RECEPTRICE
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nb1, nb2, n1, n2, i, j
    real(kind=8) :: amat(nb1, nb2), bmat(n1, n2)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ideb, ifin, ii, iideb, iifin
    integer(kind=8) :: jdeb, jfin, jj, jjdeb, jjfin
!-----------------------------------------------------------------------
    jdeb = j
    jfin = min(j+n2-1, nb2)
    if (jfin .lt. jdeb) goto 999
    jjdeb = jdeb-j+1
    jjfin = jfin-j+1
!
    ideb = i
    ifin = min(i+n1-1, nb1)
    if (ifin .lt. ideb) goto 999
    iideb = ideb-i+1
    iifin = ifin-i+1
!
    do ii = iideb, iifin
        do jj = jjdeb, jjfin
            amat(i+ii-1, j+jj-1) = amat(i+ii-1, j+jj-1)+bmat(ii, jj)
        end do
    end do
!
999 continue
end subroutine
