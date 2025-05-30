! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine xmmco3(ino, ndim, dsidep, pla, p, &
                  ffc, jac, nnol, raug, mmat)
    implicit none
#include "jeveux.h"
#include "asterfort/promat.h"
#include "asterfort/transp.h"
    integer :: ino, ndim, pla(27)
    integer :: nnol
    real(kind=8) :: mmat(216, 216), dsidep(6, 6)
    real(kind=8) :: jac, ffc(8), raug
    real(kind=8) :: p(3, 3)
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES MATRICES DE COHESION:
! --- PARTIE DEPENDANTE DE LA LOI D INTERFACE POUR CZM_LIN_MIX
!
! ----------------------------------------------------------------------
!
! IN  INO    : NUMERO LOCAL DU NOEUD
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  DSIDEP : DERIVEE DE LA CONTRAINTE COHESIVE PAR RAPPORT
!              AU MULTIPLICATEUR AUGMENTE, EN BASE LOCALE
! IN  PLA    : PLACE DES LAGRANGES DANS LA NUMEROTATION
! IN  P      : MATRICE DE CHANGEMENT DE BASE
! IN  FFC    : FONCTIONS DE FORME DE CONTACT
! IN  JAC    : PRODUIT DU JACOBIEN ET DU POIDS
! IN  NNOL   : NOMBRE DE NOEUDS PORTANT DES INCONNUES DE CONTACT
! IN  RAUG   : PARAMETRE D AUGMENTATION
! I/O MMAT   : MATRICE ELEMENTAITRE DE CONTACT/FROTTEMENT
!
!
!
!
    integer :: i, j, k, l, pli
    real(kind=8) :: dside2(3, 3), ptr(3, 3), temp(3, 3), au(3, 3)
!
! ----------------------------------------------------------------------
!
!     INITIALISATION
! on reecrit tout
! matrices A et AT
    au(:, :) = 0.d0
    dside2(:, :) = 0.d0
    temp(:, :) = 0.d0
    ptr(:, :) = 0.d0
! ensuite, on prend exemple sur XMMCO2 pour l'opération
! changement de base
    do i = 1, ndim
        do j = 1, ndim
            dside2(i, j) = dsidep(i, j)
        end do
    end do
!
! MATRICE TANGENTE EN BASE FIXE [P]T [DSIDEP] [P]
!
! ne sert plus dans la formulation avec directions
    call transp(p, 3, ndim, ndim, ptr, &
                3)
    call promat(ptr, 3, ndim, ndim, dside2, &
                3, ndim, ndim, temp)
    call promat(temp, 3, ndim, ndim, p, &
                3, ndim, ndim, au)
! on peut alors remplir la matrice C : w*/lambda et sa transposee
! on vire la boucle sur I
! elle est déjà à l'extérieur
! avec la formulation qui inclut les directions
! on prend direct la matrice locale
    pli = pla(ino)
    do l = 1, ndim
        do k = 1, ndim
            mmat(pli-1+k, pli-1+ndim+l) = mmat(pli-1+k, pli-1+ndim+l)+ffc(ino)*dside2(k, l)*jac
            mmat(pli-1+ndim+l, pli-1+k) = mmat(pli-1+ndim+l, pli-1+k)+ffc(ino)*dside2(k, l)*jac
        end do
    end do
! on remplit la matrice w*/w
    pli = pla(ino)
    do l = 1, ndim
        do k = 1, ndim
            mmat(pli-1+ndim+k, pli-1+ndim+l) = mmat(pli-1+ndim+k, pli-1+ndim+l)+raug*ffc(ino)*dsi&
                                              &de2(k, l)*jac
        end do
    end do
! on remplit la matrice lambda*/lambda :
! attention a ne pas oublier lambda*lambda*/raug
! sur la diagonale uniquement!
    pli = pla(ino)
    do l = 1, ndim
        do k = 1, ndim
            mmat(pli-1+k, pli-1+l) = mmat(pli-1+k, pli-1+l)+ffc(ino)*dside2(k, l)*jac/raug
        end do
        mmat(pli-1+l, pli-1+l) = mmat(pli-1+l, pli-1+l)-ffc(ino)*jac/raug
    end do
!
end subroutine
