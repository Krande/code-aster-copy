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
subroutine xcalf2(he, lsng, lstg, baslog, fe, &
                  dgdgl, iret)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/xbasgl.h"
#include "asterfort/xdeffe.h"
#include "asterfort/xderfe.h"
    real(kind=8) :: he, lsng, lstg, baslog(6), fe(4), dgdgl(4, 2)
    integer(kind=8) :: iret
!
!
!
!
!
!     BUT:  CALCUL DES FONCTIONS D'ENRICHISSEMENT EN UN POINT DE GAUSS
!
! IN  XYZ     : COORDONNEES DU POINT DE GAUSS CONSIDERE
! IN  HE      : VALEUR DE LA FONCTION HEAVYSIDE CSTE LE SS-ELT
! IN  LSNG    : VALEUR DE LA LEVEL SET NORMALE AU POINT DE GAUSS
! IN  LSTG    : VALEUR DE LA LEVEL SET TANGENTE AU POINT DE GAUSS
! IN  BASLOG  : BASE LOCALE AU FOND DE FISSURE AU POINT DE GAUSS
!
! OUT FE      : VALEURS DES FONCTIONS D'ENRICHISSEMENT
! OUT DGDGL   : DERIVEES DES FONCTIONS D'ENRICHISSEMENT
! OUT IRET    : CODE RETOUR VALANT 0 SI ON SE TROUVE SUR LE FOND DE
!               FISSURE (RG=0).
!               LES DERIVEES DES FONCTIONS SINGULIERES (DGDGL)
!               NE SONT ALORS PAS CALCULEES (CAR EN 1/SQRT(RG)).
!
!----------------------------------------------------------------
!
    integer(kind=8) :: i, j, k
    real(kind=8) :: p(3, 3), invp(3, 3)
    real(kind=8) :: rg, tg, dgdpo(4, 2), dgdlo(4, 2)
!
!     RECUPERATION DE LA BASE LOCALE ASSOCIEE AU PT
!     (E1=GRLT,E2=GRLN)
    call xbasgl(2, baslog, 1, p, invp)
!
!     COORDONNÉES POLAIRES DU POINT
    rg = sqrt(lsng**2+lstg**2)
!
    if (rg .gt. r8prem()) then
!       LE POINT N'EST PAS SUR LE FOND DE FISSURE
        tg = he*abs(atan2(lsng, lstg))
        iret = 1
    else
!       LE POINT EST SUR LE FOND DE FISSURE :
!       L'ANGLE N'EST PAS DÉFINI, ON LE MET À ZÉRO
!       ON NE FERA PAS LE CALCUL DES DÉRIVÉES
        tg = 0.d0
        iret = 0
    end if
!
!     FONCTIONS D'ENRICHISSEMENT DANS LA BASE POLAIRE -> FE
    call xdeffe(rg, tg, fe)
!
!     CALCUL DES DERIVEES
!     -------------------
!
!     DERIVEES DES FONCTIONS D'ENRICHISSEMENT DANS LA BASE POLAIRE
    call xderfe(rg, tg, dgdpo)
!
!     DERIVEES DES FONCTIONS D'ENRICHISSEMENT DANS LA BASE LOCALE
    do i = 1, 4
        dgdlo(i, 1) = dgdpo(i, 1)*cos(tg)-dgdpo(i, 2)*sin(tg)/rg
        dgdlo(i, 2) = dgdpo(i, 1)*sin(tg)+dgdpo(i, 2)*cos(tg)/rg
    end do
!
!     DERIVEES DES FONCTIONS D'ENRICHISSEMENT DANS LA BASE GLOBALE
    do i = 1, 4
        do j = 1, 2
            dgdgl(i, j) = 0.d0
            do k = 1, 2
                dgdgl(i, j) = dgdgl(i, j)+dgdlo(i, k)*invp(k, j)
            end do
        end do
    end do
!
end subroutine
