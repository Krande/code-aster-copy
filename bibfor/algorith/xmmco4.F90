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
subroutine xmmco4(ndim, nno, pla, nd, tau1, &
                  tau2, ffc, ddls, jac, ffp, &
                  nnol, ddlm, nnos, mmat)
    implicit none
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/transp.h"
#include "asterfort/xcalc_saut.h"
    integer :: ndim, nno, ddls, pla(27)
    integer :: nnol, ddlm, nnos
    real(kind=8) :: mmat(216, 216)
    real(kind=8) :: ffp(27), jac, ffc(8)
    real(kind=8) :: nd(3), tau1(3), tau2(3)
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES MATRICES DE COHESION, LOI CZM_LIN_MIX
! --- PARTIE INDEPENDANTE DE LA LOI D'INTERFACE
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT DE REF PARENT
! IN  PLA    : PLACE DES LAGRANGES DANS LA NUMEROTATION
! IN  ND     : DIRECTION NORMALE
! IN  TAU1   : DIRECTION TANGENTE 1
! IN  TAU2   : DIRECTION TANGENTE 2
! IN  DDLS   : NOMBRE DE DDLS DES NOEUDS SOMMET
! IN  JAC    : PRODUIT DU JACOBIEN ET DU POIDS
! IN  FFP    : FONCTIONS DE FORME DE L'ELEMENT PARENT
! IN  NNOL   : NOMBRE DE NOEUDS AVEC INCONNUES DE LAGRANGE
! IN  DDLM   : NOMBRE DE DDLS DES NOEUDS MILIEU
! IN  NNOS   : NOMBRE DE NOEUDS SOMMET
! I/O MMAT   : MATRICE ELEMENTAITRE DE COHESION
!
!
!
!
    integer :: i, j, k, l, pli, plj, jn
    real(kind=8) :: dside2(3, 3), ptr(3, 3), temp(3, 3), au(3, 3)
    real(kind=8) :: p(3, 3), coefj
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
    p(:, :) = 0.d0
    coefj = xcalc_saut(1, 0, 1)
! idem, il va falloir introduire les matrices de passage
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
! on construit la transposee de la matrice de passage
    call transp(p, 3, ndim, ndim, ptr, &
                3)
!
    do i = 1, nnol
        pli = pla(i)
        do j = 1, nno
            call indent(j, ddls, ddlm, nnos, jn)
            do l = 1, ndim
                do k = 1, ndim
! on remplit A : matrice [u*] / mu
                    mmat(pli+2*ndim-1+k, jn+ndim+l) = mmat(pli+2*ndim-1+k, jn+ndim+l)+coefj*ffc(i)&
                                                     &*p(k, l)*ffp(j)*jac
! et sa transposee
                    mmat(jn+ndim+l, pli+2*ndim-1+k) = mmat(jn+ndim+l, pli+2*ndim-1+k)+coefj*ffc(i)&
                                                     &*p(k, l)*ffp(j)*jac
                end do
            end do
        end do
    end do
!
! on remplit B : matrice w* / mu
    do i = 1, nnol
        pli = pla(i)
        do j = 1, nnol
            plj = pla(j)
            do l = 1, ndim
! on remplit B
                mmat(pli-1+ndim+l, plj+2*ndim-1+l) = mmat(pli-1+ndim+l, plj+2*ndim-1+l)-ffc(i)*ff&
                                                    &c(j)*jac
! et sa transposee
                mmat(plj+2*ndim-1+l, pli+ndim-1+l) = mmat(plj+2*ndim-1+l, pli+ndim-1+l)-ffc(i)*ff&
                                                    &c(j)*jac
!
            end do
        end do
    end do
! on a enlevé le remplissage de C
!
end subroutine
