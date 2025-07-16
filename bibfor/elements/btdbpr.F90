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
subroutine btdbpr(b, d, jacob, nbsig, nbinco, &
                  btdb)
!.======================================================================
    implicit none
!
!       BTDBPR  -- CALCUL DU PRODUIT BT*D*B DONNANT LA MATRICE
!                  DE RIGIDITE ELEMENTAIRE EN FAISANT LE PRODUIT
!                  MATRICIEL MULTIPLIE PAR LE SCALAIRE JACOBIEN*POIDS
!
!   ARGUMENT        E/S  TYPE         ROLE
!    B(NBSIG,1)     IN     R        MATRICE (B) CALCULEE AU POINT
!                                   D'INTEGRATION COURANT ET RELIANT
!                                   LES DEFORMATIONS DU PREMIER ORDRE
!                                   AUX DEPLACEMENTS
!    D(NBSIG,1)     IN     R        MATRICE DE HOOKE DANS LE REPERE
!                                   GLOBAL
!    JACOB          IN     R        PRODUIT JACOBIEN*POIDS AU POINT
!                                   D'INTEGRATION COURANT
!    NBSIG          IN     I        NOMBRE DE CONTRAINTES ASSOCIE A
!                                   L'ELEMENT
!    NBINCO         IN     I        NOMBRE D'INCONNUES SUR L'ELEMENT
!    BTDB(81,81)    OUT    R        MATRICE ELEMENTAIRE DE RIGIDITE
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
    integer(kind=8) :: nbinco, nbsig
    real(kind=8) :: b(nbsig, 1), d(nbsig, 1), jacob, btdb(81, 1)
! -----  VARIABLES LOCALES
    real(kind=8) :: tab1(10), tab2(10)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, j1, j2
    real(kind=8) :: s, zero
!-----------------------------------------------------------------------
    zero = 0.0d0
!
    do i = 1, nbinco
        do j = 1, nbsig
            tab1(j) = jacob*b(j, i)
        end do
!
        do j1 = 1, nbsig
            s = zero
            do j2 = 1, nbsig
                s = s+tab1(j2)*d(j1, j2)
            end do
            tab2(j1) = s
        end do
!
        do j1 = 1, i
            s = zero
            do j2 = 1, nbsig
                s = s+b(j2, j1)*tab2(j2)
            end do
!
            btdb(i, j1) = btdb(i, j1)+s
            btdb(j1, i) = btdb(i, j1)
!
        end do
    end do
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
