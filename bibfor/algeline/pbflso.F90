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
subroutine pbflso(umoy, rmoy, long, icoq, imod, &
                  nbm, rkip, tcoef, harm, lambda, &
                  kcalcu, passag, condit, gamma, d, &
                  ysol)
    implicit none
! COUPLAGE FLUIDELASTIQUE, CONFIGURATIONS DU TYPE "COQUE_COAX"
! RESOLUTION DU PROBLEME FLUIDE INSTATIONNAIRE : CALCUL DU VECTEUR
! SOLUTION T(UI*,VI*,PI*) LE LONG DES STRUCTURES
! APPELANT : PBFLUI
!-----------------------------------------------------------------------
!  IN : UMOY   : VITESSE DE L'ECOULEMENT MOYEN
!  IN : RMOY   : RAYON MOYEN
!  IN : LONG   : LONGUEUR DU DOMAINE DE RECOUVREMENT DES DEUX COQUES
!  IN : ICOQ   : INDICE CARACTERISANT LA COQUE SUR LAQUELLE ON TRAVAILLE
!                ICOQ=1 COQUE INTERNE  ICOQ=2 COQUE EXTERNE
!  IN : IMOD   : INDICE DU MODE CONSIDERE
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!  IN : RKIP   : ORDRE DE COQUE DU MODE CONSIDERE, PONDERE PAR LA VALEUR
!                MOYENNE DU PROFIL DE PRESSION
!  IN : TCOEF  : TABLEAU DES COEFFICIENTS DES DEFORMEES AXIALES
!  IN : HARM   : VECTEUR DE TRAVAIL
!  IN : LAMBDA : VALEURS PROPRES DE L'OPERATEUR DIFFERENTIEL
!  IN : KCALCU : MATRICE RECTANGULAIRE A COEFFICIENTS CONSTANTS
!                PERMETTANT DE CALCULER UNE SOLUTION PARTICULIERE DU
!                PROBLEME FLUIDE INSTATIONNAIRE, LORSQUE UMOY > 0
!  IN : PASSAG : MATRICE DONT LES COLONNES SONT LES VECTEURS PROPRES DE
!                L'OPERATEUR DIFFERENTIEL
!  IN : CONDIT : COEFFICIENTS DE PRECONDITIONNEMENT
!  IN : GAMMA  : COEFFICIENTS DE LA COMBINAISON LINEAIRE DONNANT LA
!                SOLUTION GENERALE DU PROBLEME FLUIDE INSTATIONNAIRE
!                (DECOMPOSITION SUR UNE FAMILLE D'EXPONENTIELLES)
!                LORSQUE UMOY > 0
!  IN : D      : COEFFICIENTS DE LA COMBINAISON LINEAIRE DONNANT LA
!                PRESSION PERTURBEE (DECOMPOSITION SUR UNE FAMILLE
!                DE FONCTIONS EXPONENTIELLES REELLES ET COMPLEXES)
!                LORSQUE UMOY = 0
! OUT : YSOL   : TABLEAU SOLUTION (VECTEUR T(UI*,VI*,PI*) TABULE EN Z)
!-----------------------------------------------------------------------
!
#include "asterfort/pbflkz.h"
    real(kind=8) :: umoy, rmoy, long
    integer(kind=8) :: icoq, imod, nbm
    real(kind=8) :: rkip, tcoef(10, nbm), harm(6)
    complex(kind=8) :: lambda(3), kcalcu(3, 4), passag(3, 3)
    real(kind=8) :: condit(3)
    complex(kind=8) :: gamma(3)
    real(kind=8) :: d(6)
    complex(kind=8) :: ysol(3, 101)
!
    real(kind=8) :: ln
    complex(kind=8) :: somm1, somm2
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: itab, k, m1, m2
    real(kind=8) :: dz, u, v, z
!-----------------------------------------------------------------------
    itab = 0
    if (icoq .eq. 2) itab = 5
    ln = tcoef(1+itab, imod)
    dz = long/100.d0
!
    if (umoy .lt. 1.d-5) then
!
        do k = 1, 101
            z = dble(k-1)*dz
            somm1 = dcmplx(0.d0, 0.d0)
            u = (ln/long)*z
            harm(1) = dble(cos(u))
            harm(2) = dble(sin(u))
            harm(3) = dble(cosh(u))
            harm(4) = dble(sinh(u))
            v = -1.d0*(rkip/rmoy)*(long-z)
            harm(5) = dble(exp(v))
            v = -1.d0*(rkip/rmoy)*z
            harm(6) = dble(exp(v))
            do m1 = 1, 6
                somm1 = somm1+d(m1)*harm(m1)
            end do
            ysol(3, k) = somm1
        end do
!
    else
!
        do k = 1, 101
            z = dble(k-1)*dz
            do m1 = 1, 3
                somm2 = dcmplx(0.d0, 0.d0)
                do m2 = 1, 3
                    somm2 = somm2+passag(m1, m2)*gamma(m2)*dcmplx(exp(lambda(m2)*(z-condit(m&
                            &2)*long)))
                end do
                ysol(m1, k) = pbflkz(m1, z, long, ln, kcalcu)+somm2
            end do
        end do
!
    end if
!
end subroutine
