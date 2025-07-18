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
subroutine dortvp(ndim, nomrc, d, modeli)
!.======================================================================
    implicit none
!
!      DORTVP --   CALCUL DES VALEURS PROPRES DE LA MATRICE
!                  HOOKE DORTH POUR S'ASSURER QUE CELLE-CI EST BIEN
!                  DEFINIE POSITIVE DANS LE CAS DE L'ORTHOTROPIE
!                  OU DE L'ISOTROPIE TRANSVERSE
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NDIM           IN    I       DIMENSION DU MODELE
!    NOMRC          IN   K16      NOM DE LA RELATION DE COMPORTEMENT
!    D(6,6)         IN    R       MATRICE DE HOOKE
!    MODELI         IN    K2      INDICATEUR DE LA MODELISATION
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterfort/infniv.h"
#include "asterfort/jacobi.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ndim
    real(kind=8) :: d(6, 6)
    character(len=2) :: modeli
    character(len=16) :: nomrc
! -----  VARIABLES LOCALES
    integer(kind=8) :: type, iordre
!
    real(kind=8) :: tr(21), tu(21), jacaux(6)
    real(kind=8) :: vecp2(4, 4), vecp3(6, 6), valp(6)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifm, ineg, k, nbvec, nitjac, niv
    integer(kind=8) :: nperm
    real(kind=8) :: tol, toldyn, un, zero
!-----------------------------------------------------------------------
    data nperm, tol, toldyn/12, 1.d-10, 1.d-2/
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     ---------------
    zero = 0.0d0
    un = 1.0d0
!
    type = 0
    iordre = 0
!
    do i = 1, 6
        valp(i) = zero
    end do
!
! --- CONSTRUCTION DU VECTEUR TR QUI CONTIENT LES VECTEURS COLONNES
! --- DE LA DEMI-MATRICE INFERIEURE DE LA MATRICE DONT ON RECHERCHE
! --- LES VALEURS PROPRES ET DU VECTEUR TU QUI CONTIENT LES VECTEURS
! --- COLONNES DE LA DEMI-MATRICE INFERIEURE DE LA MATRICE UNITE :
!     ==========================================================
!
! --- CAS 3D :
!     ------
    if (ndim .eq. 3) then
!
        nbvec = 6
!
! ---   TABLEAU TR :
!       ----------
        tr(1) = d(1, 1)
        tr(2) = d(2, 1)
        tr(3) = d(3, 1)
        tr(4) = d(4, 1)
        tr(5) = d(5, 1)
        tr(6) = d(6, 1)
!
        tr(7) = d(2, 2)
        tr(8) = d(3, 2)
        tr(9) = d(4, 2)
        tr(10) = d(5, 2)
        tr(11) = d(6, 2)
!
        tr(12) = d(3, 3)
        tr(13) = d(4, 3)
        tr(14) = d(5, 3)
        tr(15) = d(6, 3)
!
        tr(16) = d(4, 4)
        tr(17) = d(5, 4)
        tr(18) = d(6, 4)
!
        tr(19) = d(5, 5)
        tr(20) = d(6, 5)
!
        tr(21) = d(6, 6)
!
! ---   TABLEAU TU :
!       ----------
        do i = 1, 21
            tu(i) = zero
        end do
!
        k = 1
        do i = 1, 6
            tu(k) = un
            k = k+7-i
        end do
!
! --- CAS 2D :
!     ------
    else if (ndim .eq. 2) then
!
        nbvec = 4
!
! ---   TABLEAU TR :
!       ----------
        tr(1) = d(1, 1)
        tr(2) = d(2, 1)
        tr(3) = d(3, 1)
        tr(4) = d(4, 1)
!
        tr(5) = d(2, 2)
        tr(6) = d(3, 2)
        tr(7) = d(4, 2)
!
        tr(8) = d(3, 3)
        tr(9) = d(4, 3)
!
        tr(10) = d(4, 4)
!
! ---   TABLEAU TU :
!       ----------
        do i = 1, 10
            tu(i) = zero
        end do
!
        k = 1
        do i = 1, 4
            tu(k) = un
            k = k+5-i
        end do
!
    end if
!
! --- RECHERCHE DES VALEURS PROPRES DE D :
!     ==================================
    if (ndim .eq. 3) then
        call jacobi(nbvec, nperm, tol, toldyn, tr, &
                    tu, vecp3, valp, jacaux, nitjac, &
                    type, iordre)
    else if (ndim .eq. 2) then
        call jacobi(nbvec, nperm, tol, toldyn, tr, &
                    tu, vecp2, valp, jacaux, nitjac, &
                    type, iordre)
    end if
!
! --- RECUPERATION DU NIVEAU D'IMPRESSION :
!     -----------------------------------
    call infniv(ifm, niv)
!
    ineg = 0
    do i = 1, nbvec
        if (valp(i) .lt. zero) then
            ineg = ineg+1
        end if
    end do
!
    if (ineg .gt. 1) then
!
        if (modeli .eq. 'CP') then
            write (ifm, 1080)
            write (ifm, 1100)
        else if (modeli .eq. 'DP') then
            write (ifm, 1090)
            write (ifm, 1110)
        end if
!
        write (ifm, 1010)
        write (ifm, 1020) nomrc
        call utmess('A', 'MATERIAL2_20')
!
    else if (ineg .eq. 1) then
!
        if (modeli .eq. 'CP') then
            write (ifm, 1080)
            write (ifm, 1100)
        else if (modeli .eq. 'DP') then
            write (ifm, 1090)
            write (ifm, 1110)
        end if
!
        write (ifm, 1030)
        write (ifm, 1020) nomrc
        call utmess('A', 'MATERIAL2_20')
    end if
!
    if (ineg .gt. 0) then
        write (ifm, 1050)
        write (ifm, 1060)
!
        do i = 1, nbvec
            write (ifm, 1070) i, valp(i)
        end do
!
        write (ifm, 1060)
        write (ifm, 1040)
!
    end if
!
! 1000 FORMAT(7X,'VALEUR PROPRE NUMERO ',I1,' DE LA MATRICE DE HOOKE ',
!     +        E12.5)
1010 format(7x, 'LA MATRICE DE HOOKE A DES VALEURS PROPRES NEGATIVES ')
1020 format(7x, 'POUR LA RELATION DE COMPORTEMENT : ', a16,/)
1030 format(7x, 'LA MATRICE DE HOOKE A UNE VALEUR PROPRE NEGATIVE ')
1040 format(/)
1050 format(7x, 'VALEURS PROPRES DE LA MATRICE DE HOOKE : ',/)
1060 format(7x, '************************')
1070 format(7x, '!  ', i1, '  !  ', e12.5, '  !')
1080 format(7x, 'TRAITEMENT DU CAS DES CONTRAINTES PLANES :')
1100 format(7x, '---------------------------------------- ')
1090 format(7x, 'TRAITEMENT DES CAS DEFORMATIONS PLANES ET AXI :')
1110 format(7x, '--------------------------------------------- ')
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
