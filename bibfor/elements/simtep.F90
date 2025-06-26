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
subroutine simtep(fami, nno, ndim, nbsig, npg, &
                  ipoids, ivf, idfde, xyz, depl, &
                  instan, angl_naut, mater, nharm, sigma)
    implicit none
!
!      SIGVMC   -- CALCUL DES  CONTRAINTES 'VRAIES'
!                  POUR LE CALCUL DE L'ENERGIE POTENTIELLE
!                  (I.E.  1/2*SIGMA_MECA - SIGMA_THERMIQUES)
!                  AUX POINTS D'INTEGRATION POUR LES ELEMENTS
!                  ISOPARAMETRIQUES
!
!   ROUTINE IDENTIQUE A SIGVMC MAIS 1/2 POUR LE TERME DE MECANIQUE
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NNO            IN     I        NOMBRE DE NOEUDS DE L'ELEMENT
!    NDIM           IN     I        DIMENSION DE L'ELEMENT (2 OU 3)
!    NBSIG          IN     I        NOMBRE DE CONTRAINTES ASSOCIE
!                                   A L'ELEMENT
!    NPG            IN     I        NOMBRE DE POINTS D'INTEGRATION
!                                   DE L'ELEMENT
!    NI(1)          IN     R        FONCTIONS DE FORME
!    DNIDX(1)       IN     R        DERIVEES DES FONCTIONS DE FORME
!                                   / X SUR L'ELEMENT DE REFERENCE
!    DNIDY(1)       IN     R        DERIVEES DES FONCTIONS DE FORME
!                                   / Y SUR L'ELEMENT DE REFERENCE
!    DNIDZ(1)       IN     R        DERIVEES DES FONCTIONS DE FORME
!                                   / Z SUR L'ELEMENT DE REFERENCE
!    POIDS(1)       IN     R        POIDS D'INTEGRATION
!    XYZ(1)         IN     R        COORDONNEES DES CONNECTIVITES
!    DEPL(1)        IN     R        VECTEUR DES DEPLACEMENTS SUR
!                                   L'ELEMENT
!    INSTAN         IN     R        INSTANT DE CALCUL
!    ANGL_NAUT(3)   IN     R        ANGLES NAUTIQUES DEFINISSANT LE REPERE
!                                   D'ORTHOTROPIE
!    MATER          IN     I        MATERIAU
!    NHARM          IN     R        NUMERO D'HARMONIQUE
!    SIGMA(1)       OUT    R        CONTRAINTES AUX POINTS D'INTEGRATION
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "jeveux.h"
#include "asterfort/sigmmc.h"
#include "asterfort/sigtmc.h"
    integer(kind=8) :: ipoids, ivf, idfde
    character(len=*) :: fami
    real(kind=8) :: xyz(1), depl(1), angl_naut(3), sigma(1)
    real(kind=8) :: instan, nharm
! -----  VARIABLES LOCALES
    integer(kind=8) :: i, mater, nbsig, ndim, nno, npg
    character(len=16) :: k16bid
    real(kind=8) :: sigth(162), zero
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     -----------------
    zero = 0.0d0
    k16bid = ' '
!
    do i = 1, nbsig*npg
        sigma(i) = zero
    end do
!
! --- CALCUL DES CONTRAINTES MECANIQUES AUX POINTS D'INTEGRATION
!      ---------------------------------------------------------
    call sigmmc(fami, nno, ndim, nbsig, npg, &
                ipoids, ivf, idfde, xyz, depl, &
                instan, angl_naut, mater, nharm, sigma)
!
! --- CALCUL DES CONTRAINTES THERMIQUES AUX POINTS D'INTEGRATION
!      ---------------------------------------------------------
    call sigtmc(fami, ndim, nbsig, npg, &
                instan, mater, angl_naut, &
                k16bid, sigth)
!
! --- CALCUL DES CONTRAINTES TOTALES AUX POINTS D'INTEGRATION
!      ---------------------------------------------------------
    do i = 1, nbsig*npg
        sigma(i) = 0.5d0*sigma(i)-sigth(i)
    end do
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
