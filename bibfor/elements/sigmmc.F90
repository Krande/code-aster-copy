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

subroutine sigmmc(fami, nno, ndim, nbsig, npg, &
                  ipoids, ivf, idfde, xyz, depl, &
                  instan, angl_naut, mater, nharm, sigma)
!.======================================================================
    implicit none
!
!      SIGMMC   -- CALCUL DES  CONTRAINTES AUX POINTS D'INTEGRATION
!                  POUR LES ELEMENTS ISOPARAMETRIQUES
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NNO            IN     I        NOMBRE DE NOEUDS DE L'ELEMENT
!    NDIM           IN     I        DIMENSION DE L'ELEMENT (2 OU 3)
!    NBSIG          IN     I        NOMBRE DE CONTRAINTES ASSOCIE
!                                   A L'ELEMENT
!    NPG            IN     I        NOMBRE DE POINTS D'INTEGRATION
!                                   DE L'ELEMENT
!    IPOIDS         IN     I        POIDS D'INTEGRATION
!    IVF            IN     I        FONCTIONS DE FORME
!    IDFDE          IN     I        DERIVEES DES FONCTIONS DE FORME
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
#include "asterfort/bmatmc.h"
#include "asterfort/dbudef.h"
#include "asterfort/dmatmc.h"
#include "asterfort/lteatt.h"
    real(kind=8) :: xyz(1), depl(1), angl_naut(3), sigma(1)
    real(kind=8) :: instan, nharm
    character(len=*) :: fami
! -----  VARIABLES LOCALES
    integer(kind=8) :: igau, nno
    real(kind=8) :: b(486), d(36), jacgau
    character(len=2) :: k2bid
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     -----------------
!-----------------------------------------------------------------------
    integer(kind=8) :: idfde, ipoids, ivf, mater, nbinco, nbsig
    integer(kind=8) :: ndim, ndim2, npg
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    k2bid = '  '
    zero = 0.0d0
    nbinco = ndim*nno
    ndim2 = ndim
    if (lteatt('FOURIER', 'OUI')) then
        ndim2 = 2
    end if
!
    sigma(1:nbsig*npg) = zero
!
! --- CALCUL DES CONTRAINTES AUX POINTS D'INTEGRATION
! ---  BOUCLE SUR LES POINTS D'INTEGRATION
!      -----------------------------------
    do igau = 1, npg
!
!  --      CALCUL DE LA MATRICE B RELIANT LES DEFORMATIONS DU
!  --      PREMIER ORDRE AUX DEPLACEMENTS AU POINT D'INTEGRATION
!  --      COURANT : (EPS_1) = (B)*(UN)
!          ----------------------------
        call bmatmc(igau, nbsig, xyz, ipoids, ivf, &
                    idfde, nno, nharm, jacgau, b)
!
!  --      CALCUL DE LA MATRICE DE HOOKE (LE MATERIAU POUVANT
!  --      ETRE ISOTROPE, ISOTROPE-TRANSVERSE OU ORTHOTROPE)
!          -------------------------------------------------
        call dmatmc(fami, mater, instan, '+', &
                    igau, 1, angl_naut, nbsig, &
                    d)
!
!  --      CALCUL DE LA CONTRAINTE AU POINT D'INTEGRATION COURANT
!          ------------------------------------------------------
        call dbudef(depl, b, d, nbsig, nbinco, &
                    sigma(1+nbsig*(igau-1)))
    end do
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
