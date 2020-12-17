! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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

subroutine thetapdg(ndim, nno, ff, dfdi, ndimte, ithet, dtdm)
    implicit none
!
#include "jeveux.h"
!#include "asterfort/plegen.h"
!#include "asterfort/dplegen.h"
!
    integer, intent(in) :: ndim, nno
    integer, intent(in) :: ndimte, ithet
    real(kind=8), intent(in) :: ff(nno), dfdi(nno, ndim)
    real(kind=8), intent(out) :: dtdm(3, 4)
!
!.......................................................................
!
!     BUT:  CALCUL DES ELEMENTS CINEMATIQUES (MATRICES F ET E, RAYON R)
!           EN UN POINT DE GAUSS (EVENTUELLEMENT EN GRANDES TRANSFORM.)
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  FF      : FONCTIONS DE FORMES
! IN  DFDI    : DERIVEE DES FONCTIONS DE FORME
! IN  ndimte  : Nombre de fonction pour la discrétisation
!                -- ndeg si legendre
!                -- nnof si lagrange
! IN  ITHET   : zr(ithet) est tel que, pour le noeud i
!                   - zr(ithet-1+6*(i-1)+1): theta_0 évalué au noeud i
!                   - zr(ithet-1+6*(i-1)+2:ithet-1+6*(i-1)+4): la 
!                        direction de proagation t évaluée au noeud i
!                   - zr(ithet-1+6*(i-1)+5): l'abscisse curviligne s
!                        évaluée au noeud i
!                   - zr(ithet-1+6*(i-1)+6) : longueur de la fissure
! OUT DTDM    : GRADIENT DE THETA
!......................................................................
!
    integer      :: i, j, k, test
!    real(kind=8) :: gam, dgam
    real(kind=8) :: th0, s, t(3)
    real(kind=8) :: gradth0(3), grads(3), gradt(3, 3)
!
!-- Calcul de theta_0, s et t au point de Gauss
    th0 = 0.d0
    s = 0.d0
    t = 0.d0
!
!
    do i=1, nno
        th0 = th0 + ff(i)*zr(ithet-1+6*(i-1)+1)
        s = s + ff(i)*zr(ithet-1+6*(i-1)+5)
        do j=1, ndim
            t(j) = t(j) + ff(i)*zr(ithet-1+6*(i-1)+j+1)
        end do
    end do
!
!-- Calcul de grad(theta_0), grad(s) et grad(t) au point de Gauss
    gradth0 = 0.d0
    grads = 0.d0
    gradt = 0.d0
!
    do k = 1, ndim
       do i=1, nno
           gradth0(k) = gradth0(k) + dfdi(i, k)*zr(ithet-1+6*(i-1)+1)
           grads(k) = grads(k) + dfdi(i,k)*zr(ithet-1+6*(i-1)+5)
           do j=1, ndim
               gradt(j,k) = gradt(j,k) + dfdi(i, k)*zr(ithet-1+6*(i-1)+j+1)
           end do
        end do
    end do
!
! ===========================================
!                   CAS 3D: A FAIRE
! ===========================================
!-- Calcul du polynome de legendre 3D : A FAIRE
!   call plegen(zi(ideg), s, xl, gam)
!   xl est recuperable avec zr(ithet-1+6*(i-1)+6) : longueur de la fissure
!   Zi(deg) est recuperable avec ndimte qui donne le nbr de polynome
!
!   a supprimer quand programmation 3D 
    test = ndimte 
!
!-- Calcul de la derivee du polynome de legendre : A FAIRE
!    call dplegen(zi(ideg), s, xl, dgam)

!-- Il y a le cas LAGRANGE à faire aussi
! ===========================================
! ===========================================

!-- Calcul de grad(theta) "analytique" au point de Gauss
    dtdm = 0.d0
    do i=1, ndim
        do j=1, ndim
!---------- CAS 2D
            dtdm(i,j) = t(i)*gradth0(j) + th0*gradt(i,j)
!---------- CAS 3D
!           dtdm(i,j) = t(i)*(dgam*th0*grads(j) + gam*gradth0(j)) + gam*th0*gradt(i, j)
        enddo
    enddo
!
!-- Stockage de theta dans la quatrième colonne de dtdm
    do i=1, ndim
!------ CAS 2D
        dtdm(i,4) = th0*t(i)
!------ CAS 3D
!        dtdm(i,4) = gam*th0*t(i)
    enddo
!
end subroutine
