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
subroutine dfdm3j(nno, ipg, idfde, coor, jac)
    implicit none
#include "jeveux.h"
    integer(kind=8) :: ipg, idfde, nno
    real(kind=8) :: coor(1), jac
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DU JACOBIEN (AVEC SIGNE)
!               POUR LES ELEMENTS 3D
!
!    - ARGUMENTS:
!        DONNEES:     NNO           -->  NOMBRE DE NOEUDS
!              DFDRDE,DFRDN,DFRDK   -->  DERIVEES FONCTIONS DE FORME
!                     COOR          -->  COORDONNEES DES NOEUDS
!
!        RESULTATS:  JAC           <--  JACOBIEN AU POINT DE GAUSS
! ......................................................................
!
    integer(kind=8) :: i, j, ii, k
    real(kind=8) :: g(3, 3)
    real(kind=8) :: de, dn, dk, j11, j21, j31
!
!
    do i = 1, 3
        do j = 1, 3
            g(i, j) = 0.d0
        end do
    end do
!
    do i = 1, nno
        k = 3*nno*(ipg-1)
        ii = 3*(i-1)
        de = zr(idfde-1+k+ii+1)
        dn = zr(idfde-1+k+ii+2)
        dk = zr(idfde-1+k+ii+3)
        do j = 1, 3
            g(1, j) = g(1, j)+coor(ii+j)*de
            g(2, j) = g(2, j)+coor(ii+j)*dn
            g(3, j) = g(3, j)+coor(ii+j)*dk
        end do
    end do
!
    j11 = g(2, 2)*g(3, 3)-g(2, 3)*g(3, 2)
    j21 = g(3, 1)*g(2, 3)-g(2, 1)*g(3, 3)
    j31 = g(2, 1)*g(3, 2)-g(3, 1)*g(2, 2)
!
    jac = g(1, 1)*j11+g(1, 2)*j21+g(1, 3)*j31
end subroutine
