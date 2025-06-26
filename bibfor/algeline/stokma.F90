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
subroutine stokma(amat, nlig, ncol, nmat, amatst)
!
! FONCTION: COPIE LA MATRICE AMAT(I,J), A NLIG LIGNES ET NCOL COLONNES,
!           DANS LA MATRICE AMATST(I,J,NMAT)
!
!     IN  : AMAT      : MATRICE DIMENSIONNEE NLIG,NCOL
!           NLIG      : ENTIER
!           NCOL      : ENTIER
!           NMAT      : ENTIER
!
!     OUT : AMATST    : TABLEAU DIMENSIONNE 9,6,6
! ------------------------------------------------------------------
    implicit none
    integer(kind=8) :: ncol, nlig, nmat
    real(kind=8) :: amat(nlig, ncol), amatst(9, 6, 6)
    integer(kind=8) :: i, j
!-----------------------------------------------------------------------
!
    do j = 1, ncol
        do i = 1, nlig
            amatst(i, j, nmat) = amat(i, j)
        end do
    end do
end subroutine
