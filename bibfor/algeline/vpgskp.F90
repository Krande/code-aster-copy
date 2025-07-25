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
subroutine vpgskp(nbeq, nconv, vect, alpha, lmatb, &
                  typeps, vaux, ddlexc, delta)
!
!     SUBROUTINE ASTER ORTHONORMALISANT LES NCONV VECTEUR VECT(:,J)
! VIA L'ALGORITHME DE TYPE GRAM-SCHMIDT ITERATIF (IGS) SUIVANT LA
! VERSION DE KAHAN-PARLETT.
!---------------------------------------------------------------------
! L'ORTHONORMALISATION EST REALISEE AU SENS DU PRODUIT SCALAIRE:
!                             (MATRICE (B) DOIT ETRE SYMETRIQUE + ...)
!   SI TYPEPS=0 -> L2         (QUELCONQUE -> PRODUIT SCALAIRE EUCLIDIEN)
!   SI TYPEPS=1 -> LMATB*SIGN (INDEFINIE -> PSEUDO-PRODUIT SCALAIRE)
!   SI TYPEPS=2 -> LMATB  (SEMI DEFINIE POSITIVE -> B-PRODUIT SCALAIRE)
!   -------------------------------------------------------------------
!     PARTANT DE LA MATRICE DE VECTEUR (V1,V2..VI..VN)
!     ON CHERCHE A OBTENIR LA MATRICE ORTHOGONALE (Q1,Q2..QN)
!
!     (0) NORMALISATION DE V1: Q1 <- V1/||V1||
!     (I) BOUCLE I = 2..N
!     (IJ)  BOUCLE J = 1..I-1
!             CALCUL DE (VJ,VI)
!             CALCUL DE ||VI||
!             CALCUL DE VI+ <- VI - (VJ,VI)VJ
!             CALCUL DE ||VI+||
!             SI ||VI+|| >= ALPHA * ||VI|| ALORS (TEST 1)
!               QI <- VI+
!             SINON
!               CALCUL DE (VJ,VI+)
!               CALCUL DE VI++ <- VI+ - (VJ,VI+)VJ
!               CALCUL DE ||VI++||
!               SI ||VI++|| >= ALPHA * ||VI+|| ALORS (TEST 2)
!                 QI <- VI++
!               SINON
!                 QI <- 0
!               FIN TEST 2
!             FIN DE TEST 1
!           FIN DE BOUCLE J
!         FIN DE BOUCLE I
!   -------------------------------------------------------------------
!     SUBROUTINES APPELLEES:
!
!       MRMULT -> (SUBROUTINE ASTER)
!         PRODUIT MATRICE-VECTEUR.
!       R8MIEM -> (FONCTION ASTER)
!         LA VALEUR MINIMUM EVITANT L'OVERFLOW LORSQU'ON L'INVERSE.
!
!     FONCTIONS INTRINSEQUES:
!
!       MAX, SIGN, ABS, SQRT.
!   --------------------------------------------------------------------
!     PARAMETRES D'APPELS:
!
! IN  NBEQ   : IS : DIMENSION DES VECTEURS.
! IN  NCONV  : IS : NOMBRE DE VALEURS PROPRES CONVERGEES.
! IN  ALPHA  : R8 : PARAMETRE DE L'ALGORITHME DE KAHAN-PARLETT.
! IN  LMATB  : IS : DESCRIPTEUR MATRICE DE PRODUIT SCALAIRE.
! IN  TYPEPS : IS : TYPE DE PRODUIT SCALAIRE.
! IN  DDLEXC : IS : DDLEXC(1..NBEQ) VECTEUR POSITION DES DDLS BLOQUES.
!
! IN/OUT VECT  : R8 : VECT(1..NBEQ,1..NCONV) MATRICE DES VECTEURS A
!                     ORTHONORMALISER (IN),
!                     MATRICE ORTHOGONALE RESULTAT (OUT).
! IN/OUT VAUX  : R8 : VAUX(1..NBEQ) VECTEUR DE TRAVAIL.
! IN/OUT DELTA : R8 : DELTA(1..NCONV) VECTEUR DE STOCKAGE DES SIGN(PS).
!
! ASTER INFORMATION
! 11/01/2000 TOILETTAGE DU FORTRAN SUIVANT LES REGLES ASTER,
!            REMPLACEMENT DE DLAMCH PAR R8MIEM,
!            REMPLACEMENT DE DABS, DSQRT, DSIGN.
!--------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterc/r8miem.h"
#include "asterfort/mrmult.h"
    integer(kind=8) :: nbeq, nconv, lmatb, typeps, ddlexc(nbeq)
    real(kind=8) :: vect(nbeq, nconv), alpha, vaux(nbeq), delta(nconv)
!
!----------------------------------------------------------------
! DECLARATION VARIABLES LOCALES
!
    integer(kind=8) :: i, j, k, step
    real(kind=8) :: eps, raux, rauold, delta1
!----------------------------------------------------------------
! INITIALISATION DU PLUS PETIT REEL*8 EVITANT L'OVERFLOW
    eps = r8miem()**(2.0d+0/3.0d+0)
!
! NORMALISATION DU PREMIER VECTEUR PROPRE
    if (typeps .ne. 0) then
        call mrmult('ZERO', lmatb, vect(1, 1), vaux, 1, &
                    .false._1)
    else
        do i = 1, nbeq
            vaux(i) = vect(i, 1)
        end do
    end if
    raux = 0.d0
    do i = 1, nbeq
        raux = raux+vect(i, 1)*vaux(i)*ddlexc(i)
    end do
    if (typeps .eq. 1) then
        delta(1) = sign(1.d0, raux)
    else
        delta(1) = 1.d0
    end if
    raux = delta(1)/max(eps, sqrt(abs(raux)))
    do i = 1, nbeq
        vect(i, 1) = vect(i, 1)*raux*ddlexc(i)
    end do
!
! BOUCLE 1 SUR LES VECTEURS PROPRES
    do i = 2, nconv
!
! BOUCLE 2 SUR LES VECTEURS PROPRES
        do j = 1, i-1
!
! CALCUL (VJ,VI) ET ||VI|| (STEP 1)/ (VJ,VI+) (STEP 2 SI NECESSAIRE)
            step = 0
20          continue
            step = step+1
!
            if (typeps .ne. 0) then
                call mrmult('ZERO', lmatb, vect(1, i), vaux, 1, &
                            .false._1)
            else
                do k = 1, nbeq
                    vaux(k) = vect(k, i)
                end do
            end if
            raux = 0.d0
            do k = 1, nbeq
                raux = raux+vect(k, j)*vaux(k)*ddlexc(k)
            end do
            if (step .eq. 1) then
                rauold = 0.d0
                do k = 1, nbeq
                    rauold = rauold+vect(k, i)*vaux(k)*ddlexc(k)
                end do
            end if
!
! CALCUL VI+ <- VI - (VJ,VI)VJ (STEP 1)
! CALCUL VI++ <- VI+ - (VJ,VI+)VJ (STEP 2)
            delta1 = raux*delta(j)
            do k = 1, nbeq
                vect(k, i) = (vect(k, i)-delta1*vect(k, j))*ddlexc(k)
            end do
!
! CALCUL DE ||VI+|| (STEP 1) ET ||VI++|| (STEP 2)
            if (typeps .ne. 0) then
                call mrmult('ZERO', lmatb, vect(1, i), vaux, 1, &
                            .false._1)
            else if (typeps .eq. 0) then
                do k = 1, nbeq
                    vaux(k) = vect(k, i)
                end do
            end if
            raux = 0.d0
            do k = 1, nbeq
                raux = raux+vect(k, i)*vaux(k)*ddlexc(k)
            end do
!
! PREMIER TEST
            if ((sqrt(abs(raux)) .gt. (alpha*sqrt(abs(rauold))-eps)) .and. (step .le. 2)) then
                goto 60
            else if (step .eq. 1) then
                rauold = raux
                goto 20
            else if (step .eq. 2) then
                do k = 1, nbeq
                    vect(k, i) = 0.d0
                end do
            end if
!
60          continue
            if (typeps .eq. 1) then
                delta(i) = sign(1.d0, raux)
            else
                delta(i) = 1.d0
            end if
            raux = delta(i)/max(eps, sqrt(abs(raux)))
            do k = 1, nbeq
                vect(k, i) = vect(k, i)*ddlexc(k)*raux
            end do
        end do
!
    end do
!
! FIN ROUTINE VPGSKP
end subroutine
