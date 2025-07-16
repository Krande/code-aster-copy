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

subroutine btsig(lonlig, loncol, jacgau, bmat, sigma, bsigma)
    implicit none
    integer(kind=8), intent(in) :: loncol, lonlig
    real(kind=8), intent(in) :: jacgau, bmat(loncol, 81), sigma(1)
    real(kind=8), intent(out) :: bsigma(1)
!
    integer(kind=8) :: i, j
    real(kind=8) :: valbsi
!-----------------------------------------------------------------------
! --- CALCUL DU PRODUIT (BT)*(SIGMA) ,
! --- AVEC LES NOTATIONS DE LA ROUTINE , CA DONNE :
! ---       (BSIGMA) = (BMAT)*(SIGMA)*JACGAU
!     ------------------------------------------------------------------
!     IN  LONLIG  : LONGUEUR D'UNE LIGNE DE (B), SOIT NBNO*NBDDL
!     IN  LONCOL  : LONGUEUR D'UNE COLONNE DE (B), SOIT NBSIG
!     IN  JACGAU  : PRODUIT DU JACOBIEN PAR LE POIDS AU POINT
!                   D'INTEGRATION COURANT
!     IN  BMAT    : MATRICE (B) AU POINT D'INTEGRATION COURANT
!     IN  SIGMA   : VECTEUR DES CONTRAINTES AU POINT D'INTEGRATION
!                   COURANT
!     OUT BSIGMA  : VECTEUR (BT)*(SIGMA)*JACGAU
!     ------------------------------------------------------------------
!
    do i = 1, lonlig
        valbsi = 0.0d0
        do j = 1, loncol
            valbsi = valbsi+bmat(j, i)*sigma(j)
        end do
!
        bsigma(i) = bsigma(i)+valbsi*jacgau
    end do
!
end subroutine
