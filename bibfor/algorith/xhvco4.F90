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

subroutine xhvco4(ino, ndim, sigma, lamb, pla, &
                  jac, ffc, p, raug, vect)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/prmave.h"
#include "asterfort/transp.h"
    integer :: ino
    integer :: ndim
    integer :: pla(27)
    real(kind=8) :: vect(560), sigma(6)
    real(kind=8) :: jac, lamb(3), raug
    real(kind=8) :: ffc(16), p(3, 3)
!
! ======================================================================
! person_in_charge: daniele.colombo at ifpen.fr
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES SECONDS MEMBRES DE COHESION
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT DE REF PARENT
! IN  SIGMA  : VECTEUR CONTRAINTE EN REPERE LOCAL
! IN  NFH    : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  DDLS   : NOMBRE DE DDL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
! IN  JAC    : PRODUIT DU JACOBIEN ET DU POIDS
! IN  FFP    : FONCTIONS DE FORME DE L'ELEMENT PARENT
! IN  SINGU  : 1 SI ELEMENT SINGULIER, 0 SINON
! IN  RR     : DISTANCE AU FOND DE FISSURE
! I/O VTMP   : VECTEUR ELEMENTAIRE DE CONTACT/FROTTEMENT
!
!
!
!
    integer :: k, pli, ier
    real(kind=8) :: ptr(3, 3), sigglo(3), ffi
!
! ---------------------------------------------------------------------
!
! DIRECTION DU SAUT DE DEPLACEMENT TANGENT
!
! on reecrit tout
! on suppose qu on a acces à LAMB(3):
!    valeur de lambda au point de gauss
! et qu on a accès à WSAUT(3)
!    valeur de w au point de Gauss
! remplissage L1l
    ptr(:, :) = 0.d0
    sigglo(:) = 0.d0
! remplissage L2w
! on transpose la matrice
    call transp(p, 3, ndim, ndim, ptr, &
                3)
! produit matrice vecteur
    call prmave(0, ptr, 3, ndim, ndim, &
                sigma, ndim, sigglo, ndim, ier)
! attention, cette fois on ne doit pas virer
! la boucle sur les noeuds,
! mais on remplit tout avec la contribution de J
! vecteur w*
    pli = pla(ino)
    ffi = ffc(ino)
!
    do k = 1, ndim
        vect(pli-1+3+ndim+k) = vect(pli-1+3+ndim+k)+sigma(k)*ffi*jac
    end do
! vecteur lambda*
! attention a ne pas oublier le terme lambda/r
    pli = pla(ino)
    do k = 1, ndim
        vect(pli-1+3+k) = vect(pli-1+3+k)+(sigma(k)-lamb(k))*ffi*jac/raug
    end do
!
end subroutine
