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
subroutine xmvco5(ndim, nno, nnol, pla, nd, &
                  tau1, tau2, mu, ddls, jac, &
                  ffc, ffp, nnos, ddlm, wsaut, &
                  saut, vtmp)
!
!
    implicit none
#include "asterfort/indent.h"
#include "asterfort/prmave.h"
#include "asterfort/transp.h"
#include "asterfort/xcalc_saut.h"
    integer :: ndim, nno, nnol, ddlm
    integer :: ddls, pla(27)
    integer :: nnos
    real(kind=8) :: vtmp(400)
    real(kind=8) :: ffp(27), jac
!
    real(kind=8) :: mu(3), wsaut(3), saut(3)
    real(kind=8) :: ffc(8)
    real(kind=8) :: nd(3), tau1(3), tau2(3)
!
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES SECONDS MEMBRES DE COHESION, LOI CZM_LIN_MIX
! --- PARTIE INDEPENDANTE DE LA LOI D'INTERFACE
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT DE REF PARENT
! IN  NNOL   : NOMBRE DE NOEUDS LAGRANGE
! IN  PLA    : PLACE DES INCONNUS DE CONTACT DANS LA NUMEROTATION
! IN  ND     : DIRECTION NORMALE
! IN  TAU1   : DIRECTION TANGENTE 1
! IN  TAU2   : DIRECTION TANGENTE 2
! IN  MU     : VECTEUR DE MEME NOM (FORCES D'INTERFACE MOYENNES)
! IN  DDLS   : NOMBRE DE DDLS DES NOEUDS SOMMET
! IN  JAC    : PRODUIT DU JACOBIEN ET DU POIDS
! IN  FFC    : FONCTIONS DE FORME DE CONTACT
! IN  FFP    : FONCTIONS DE FORME DE L'ELEMENT PARENT
! IN  NNOS   : NOMBRE DE NOEUDS SOMMET
! IN  DDLM   : NOMBRE DE DDLS DES NOEUDS MILIEU
! IN  SAUT   : SAUT DE DEPLACEMENT (BASE FIXE)
! IN  WSAUT  : SAUT DE DEPLACEMENT MOYEN W (BASE LOCAL)
! I/O VTMP   : VECTEUR ELEMENTAIRE DE COHESION
!
!
!
!
    integer :: i, j, pli, iin, ier
    real(kind=8) :: p(3, 3), ptr(3, 3), mug(3), am(3), coefi
!
! ---------------------------------------------------------------------
! on va commencer par construire une matrice de passage
    p(:, :) = 0.d0
    ptr(:, :) = 0.d0
    mug(:) = 0.d0
    am(:) = 0.d0
    coefi = xcalc_saut(1, 0, 1)
!
    do i = 1, ndim
        p(1, i) = nd(i)
    end do
    do i = 1, ndim
        p(2, i) = tau1(i)
    end do
    if (ndim .eq. 3) then
        do i = 1, ndim
            p(3, i) = tau2(i)
        end do
    end if
! calcul du saut de deplacement am en base locale
    call prmave(0, p, 3, ndim, ndim, &
                saut, ndim, am, ndim, ier)
! calcul de la contrainte MUG en base globale
    call transp(p, 3, ndim, ndim, ptr, &
                3)
    call prmave(0, ptr, 3, ndim, ndim, &
                mu, ndim, mug, ndim, ier)
! on reecrit tout
! on suppose qu on a acces à LAMB(3):
!    valeur de lambda au point de gauss
! et qu on a accès à WSAUT(3)
!    valeur de w au point de Gauss
! remplissage L1mu
    do i = 1, nno
        call indent(i, ddls, ddlm, nnos, iin)
        do j = 1, ndim
            vtmp(iin+ndim+j) = vtmp(iin+ndim+j)+coefi*mug(j)*ffp(i)*jac
        end do
    end do
! remplissage L1u
    do i = 1, nnol
        pli = pla(i)
        do j = 1, ndim
            vtmp(pli+2*ndim-1+j) = vtmp(pli+2*ndim-1+j)-am(j)*ffc(i)*jac
        end do
    end do
! remplissage L1w
    do i = 1, nnol
        pli = pla(i)
        do j = 1, ndim
            vtmp(pli+2*ndim-1+j) = vtmp(pli+2*ndim-1+j)-wsaut(j)*ffc(i)*jac
        end do
    end do
! remplissage L2mu
    do i = 1, nnol
        pli = pla(i)
        do j = 1, ndim
            vtmp(pli-1+ndim+j) = vtmp(pli-1+ndim+j)-mu(j)*ffc(i)*jac
        end do
    end do
! le remplissage de L2w a été enlevé
!
end subroutine
