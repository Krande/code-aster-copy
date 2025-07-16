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

subroutine thetapdg(ndim, nno, discr, ff, dfdi, ideg, ilag, ithet, dtdm)
    implicit none
!
#include "asterfort/assert.h"
#include "jeveux.h"
#include "asterfort/plegen.h"
#include "asterfort/dplegen.h"
!
    integer(kind=8), intent(in) :: ndim, nno
    integer(kind=8), intent(in) :: ithet, ideg, ilag
    real(kind=8), intent(in)  :: ff(nno), dfdi(nno, ndim)
    real(kind=8), intent(out) :: dtdm(3, 4)
    character(len=8) :: discr
!
!.......................................................................
!
!     BUT:  CALCUL DES ELEMENTS CINEMATIQUES (MATRICES F ET E, RAYON R)
!           EN UN POINT DE GAUSS (EVENTUELLEMENT EN GRANDES TRANSFORM.)
!
! IN  NDIM      : DIMENSION DE L'ESPACE
! IN  NNO       : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  DISCR     : CHOIX DE LA DISCRETISATION : LAGRANGE OU LEGENDRE OU 2D
! IN  FF        : FONCTIONS DE FORMES
! IN  DFDI      : DERIVEE DES FONCTIONS DE FORME
! IN  IDEG      : Degré du polynome pour la discrétisation legendre k
! IN  ilag      : Abs_cur s0, s1 et s2 pour la fonction de lagrange k
! IN  ITHET     : zr(ithet) est tel que, pour le noeud i
!                   - zr(ithet-1+6*(i-1)+1): theta_0 évalué au noeud i
!                   - zr(ithet-1+6*(i-1)+2:ithet-1+6*(i-1)+4): la
!                        direction de propagation t évaluée au noeud i
!                   - zr(ithet-1+6*(i-1)+5): l'abscisse curviligne s
!                        évaluée au noeud i
!                   - zr(ithet-1+6*(i-1)+6) : longueur de la fissure
! OUT DTDM      : GRADIENT DE THETA + THETA (4ème colonne)
!......................................................................
!
    integer(kind=8)      :: i, j, k
    real(kind=8) :: th0, s, t(3), xl
    real(kind=8) :: sno, s0, s1, s2
    real(kind=8) :: gam, dgam, eval
    real(kind=8) :: gradth0(3), grads(3), gradt(3, 3)
!
!-- Initilalisation des paramètres
    th0 = 0.d0
    s = 0.d0
    t = 0.d0
    gradth0 = 0.d0
    grads = 0.d0
    gradt = 0.d0
    gam = 0.d0
    dgam = 0.d0
    dtdm = 0.d0
!
!-- Calcul de theta_0, s et t au point de Gauss
    do i = 1, nno
        th0 = th0+ff(i)*zr(ithet-1+6*(i-1)+1)
        s = s+ff(i)*zr(ithet-1+6*(i-1)+5)
        do j = 1, ndim
            t(j) = t(j)+ff(i)*zr(ithet-1+6*(i-1)+j+1)
        end do
    end do
!
!-- Calcul de grad(theta_0), grad(s) et grad(t) au point de Gauss
    do k = 1, ndim
        do i = 1, nno
            gradth0(k) = gradth0(k)+dfdi(i, k)*zr(ithet-1+6*(i-1)+1)
            grads(k) = grads(k)+dfdi(i, k)*zr(ithet-1+6*(i-1)+5)
            do j = 1, ndim
                gradt(j, k) = gradt(j, k)+dfdi(i, k)*zr(ithet-1+6*(i-1)+j+1)
            end do
        end do
    end do
!
! ===========================================
!          CAS 3D : discrétisation
! ===========================================
!
    if (discr == "LEGENDRE") then
!
!       Longueur de la fissure
        xl = zr(ithet-1+6*(1-1)+6)
!
!       Polynome de Legendre
        call plegen(zi(ideg), s, xl, gam)
!
!       Derivee du polynome de legendre
        call dplegen(zi(ideg), s, xl, dgam)
!
!------ Calcul de grad(theta) "analytique" au point de Gauss
!-------Stockage de theta dans la 4ème colonne
        do i = 1, ndim
            dtdm(i, 4) = gam*th0*t(i)
            do j = 1, ndim
                dtdm(i, j) = t(i)*(dgam*th0*grads(j)+gam*gradth0(j))+gam*th0*gradt(i, j)
            end do
        end do
!
    else if (discr == "LINEAIRE") then
!
!       Longueur de la fissure
        xl = zr(ithet-1+6*(1-1)+6)
!
        do i = 1, nno
!
!---------- Récupération des asbcisses curvilignes pour les noeuds du fond
!---------- restreint a la foncction de forme courante
            s0 = zr(ilag)
            s1 = zr(ilag+1)
            s2 = zr(ilag+2)
            sno = zr(ithet-1+6*(i-1)+5)
!
!---------- Recherche de la valeur de l'abscisse curviligne normalisée
!---------- pour le projeté d'un noeud du maillage
            if (sno .ge. s0 .and. sno .le. s1 .and. s0 .ne. s1) then
                eval = (sno-s0)/(s1-s0)
            elseif (sno .ge. s1 .and. sno .le. s2) then
                eval = 1.d0-((sno-s1)/(s2-s1))
!---------- Cas fond fermé
            elseif (s0 .gt. s1 .and. sno .ge. s0) then
                eval = (sno-s0)/(xl-s0)
            else
                eval = 0.d0
            end if
!
!---------  Calcul de grad(theta) au point de Gauss
!---------- Stockage de theta dans la 4ème colonne
            do j = 1, ndim
                dtdm(j, 4) = dtdm(j, 4)+ff(i)*(eval*zr(ithet-1+6*(i-1)+1)* &
                                               zr(ithet-1+6*(i-1)+j+1))
                do k = 1, ndim
                    dtdm(j, k) = dtdm(j, k)+dfdi(i, k)*(eval*zr(ithet-1+6*(i-1)+1)* &
                                                        zr(ithet-1+6*(i-1)+j+1))
                end do
            end do
        end do
!
    else if (discr == "2D") then
!
!------ Calcul de grad(theta) "analytique" au point de Gauss
!-------Stockage de theta dans la 4ème colonne
        do i = 1, ndim
            dtdm(i, 4) = th0*t(i)
            do j = 1, ndim
                dtdm(i, j) = t(i)*gradth0(j)+th0*gradt(i, j)
            end do
        end do
!
    else
!
        ASSERT(ASTER_FALSE)
!
    end if

!
end subroutine
