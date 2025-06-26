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

subroutine sigtmc(fami, ndim, nbsig, npg, &
                  instan, mater, angl_naut, &
                  option, sigma)
!.======================================================================
    implicit none
!
!      SIGTMC   -- CALCUL DES  CONTRAINTES THERMIQUES/HYDRIQUE OU DE
!                  SECHAGE AUX POINTS D'INTEGRATION
!                  POUR LES ELEMENTS ISOPARAMETRIQUES
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NDIM           IN     I        DIMENSION DE L'ELEMENT (2 OU 3)
!    NBSIG          IN     I        NOMBRE DE CONTRAINTES ASSOCIE
!                                   A L'ELEMENT
!    NPG            IN     I        NOMBRE DE POINTS D'INTEGRATION
!                                   DE L'ELEMENT
!    INSTAN         IN     R        INSTANT DE CALCUL (0 PAR DEFAUT)
!    MATER          IN     I        MATERIAU
!    ANGL_NAUT(3)   IN     R        ANGLES NAUTIQUES DEFINISSANT LE REPERE
!                                   D'ORTHOTROPIE
!    OPTION         IN     K16      OPTION DE CALCUL
!    SIGMA(1)       OUT    R        CONTRAINTES THERMIQUES
!                                   AUX POINTS D'INTEGRATION
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterc/r8miem.h"
#include "asterfort/dmatmc.h"
#include "asterfort/epstmc.h"
#include "asterfort/lteatt.h"
    character(len=16) :: option
    real(kind=8) :: angl_naut(3), sigma(1)
    real(kind=8) :: instan
    character(len=*) :: fami
! -----  VARIABLES LOCALES
    real(kind=8) :: d(36), epsth(6)
    integer(kind=8) :: iepsv
    character(len=2) :: k2bid
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     -----------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, igau, j, mater, nbsig, ndim
    integer(kind=8) :: ndim2, npg
    real(kind=8) :: zero
!-----------------------------------------------------------------------
    k2bid = '  '
    zero = 0.0d0
    ndim2 = ndim
    if (lteatt('FOURIER', 'OUI')) then
        ndim2 = 2
    end if
!
    sigma(1:nbsig*npg) = zero

!
! --- CALCUL DES CONTRAINTES D'ORIGINE THERMIQUE/HYDRIQUE/SECHAGE
! ---  BOUCLE SUR LES POINTS D'INTEGRATION
!      -----------------------------------
    do igau = 1, npg
!
!       REM : HYDRATATION ET SECHAGE ACTIVES POUR CALCUL DU SECOND
!        MEMBRE CHAR_MECA_* OU CALCUL DES CONTRAINTES VRAIES (SIGVMC.F)
!
!
!  --      CALCUL DES DEFORMATIONS THERMIQUES/HYDRIQUE/DE SECHAGE
!  --      AU POINT D'INTEGRATION COURANT
!          ------------------------------
        call epstmc(fami, ndim, instan, '+', igau, &
                    1, angl_naut, mater, option, &
                    epsth)
!
! TEST DE LA NULLITE DES DEFORMATIONS DUES AUX VARIABLES DE COMMANDE
        iepsv = 0
        do i = 1, 6
            if (abs(epsth(i)) .gt. r8miem()) iepsv = 1
        end do
! TOUTES DES COMPOSANTES SONT NULLES. ON EVITE LE CALCUL DE D ET SIGMA
        if (iepsv .ne. 0) then
!
!         PASSAGE DES COMPOSANTES DE CISAILLEMENTS EN CONFORMITE
!         ( DMATMC RETOURNE UNE MATRICE DE HOOKE EN SUPPOSANT
!           UN DEUX SUR LES DEFORMATIONS DE CISAILLEMENT )
!
            do i = 4, 2*ndim
                epsth(i) = 2.d0*epsth(i)
            end do
!  --      CALCUL DE LA MATRICE DE HOOKE (LE MATERIAU POUVANT
!  --      ETRE ISOTROPE, ISOTROPE-TRANSVERSE OU ORTHOTROPE)
!          -------------------------------------------------
            call dmatmc(fami, mater, instan, '+', igau, &
                        1, angl_naut, nbsig, d)
!
!  --      CONTRAINTES THERMIQUES/HYDRIQUE/DE SECHAGE AU POINT
!  --      D'INTEGRATION COURANT
!          ------------------------------------------------------
            do i = 1, nbsig
                do j = 1, nbsig
                    sigma(i+nbsig*(igau-1)) = sigma(i+nbsig*(igau-1))+d(j+(i-1)*nbsig)*epsth(j)
                end do
            end do
        end if
    end do
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
