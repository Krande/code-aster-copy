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
subroutine xsifl2(basloc, coeff, coeff3, ddld, ddlm,&
                  ddls, dfdi, ff, idepl, igthet,&
                  ithet, jac, ndim, nnop, nnos,&
                  tau1, tau2, nd, xg)
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/indent.h"
#include "asterfort/normev.h"
#include "asterfort/prmave.h"
#include "asterfort/provec.h"
#include "asterfort/transp.h"
#include "blas/ddot.h"
!
! Calcul de facteurs d'intensité des contraintes équivalents
!   avec XFEM + éléments cohésifs (option CALC_K_G_COHE)
!
! In basloc => base covariante liée au front,
!              aux noeuds du maillage parent
! In coeff, coeff3 => coefficients pour application de la
!                     formule d'Irwin
! In ddld => nombre de ddl de déplacement par noeud
! In ddlm => nombre de ddl par noeud milieu
! In ddls => nombre de ddl par noeud sommet
! In dfdi => dérivées des fonctions de forme
! In ff   => fonctions de forme
! In idepl => champ de déplacement
! Out igthet => champ resultat (KI, KII etc...)
! In  ithet => champ theta (extension virtuelle de fissure)
! In  jac => produit du jacobien et du poids
! In ndim => dimension
! In nnop => nombre de noeuds
! In nnos => nombre de noeuds sommet
! In nd => normale à la facette
! In tau1 => tangente 1
! In xg => coordonnées du point de Gauss
!
!
!
    integer :: ndim, nnop
    real(kind=8) :: am(3), basloc(9*nnop), coeff, coeff3
    integer :: ddld, ddlm, ddls
    real(kind=8) :: dfdi(nnop, ndim)
    real(kind=8) :: e1(3), e2(3), e3(3), ff(27), g, g1, g2, g3
    real(kind=8) :: grde1, grde2, grde3, grdep(3, 3)
    real(kind=8) :: gs2, gs3, cmp_hp
    integer :: i, idepl, ier, igthet, ii, ino, ithet
    integer :: j, l
    real(kind=8) :: jac, k1, k2, k3, lamb(3), lamb1, lamb2
    real(kind=8) :: lamb3, lambl(3), norme_theta
    integer :: nnos
    real(kind=8) :: norme, pm(3, 3), ptr(3, 3), theta(3)
    real(kind=8) :: tau1(3), tau2(3), nd(3), temp(3), xg(3)
    real(kind=8) :: ptp(3), vec(3), sens
    blas_int :: b_incx, b_incy, b_n
!
!     BASE LOCALE ET LEVEL SETS AU POINT DE GAUSS
!     DIMENSIONNEMENT A 3 ET NON NDIM POUR POUVOIR UTILISER NORMEV.F
    pm(:, :) = 0.d0
    do i = 1, ndim
        pm(1, i) = nd(i)
    end do
    do i = 1, ndim
        pm(2, i) = tau1(i)
    end do
    if (ndim .eq. 3) then
        do i = 1, ndim
            pm(3, i) = tau2(i)
        end do
    end if
    call transp(pm, ndim, ndim, ndim, ptr,&
                ndim)
    e1(:) = 0.d0
    e2(:) = 0.d0
    ptp(:) = 0.d0
    vec(:) = 0.d0
    do ino = 1, nnop
        do i = 1, ndim
            ptp(i) = ptp(i)+basloc(3*ndim*(ino-1)+i)*ff(ino)
            e2(i) = e2(i)+basloc(3*ndim*(ino-1)+i+ndim)*ff(ino)
            e1(i) = e1(i)+basloc(3*ndim*(ino-1)+i+2*ndim)*ff(ino)
        end do
    end do
!
! E1 == NORMALE (N)
! E2 == TANGENTE DIRECTION DE FISSURATION (M)
!
!     -----------------------------------
!     2) CALCUL DE THETA ET DE DIV(THETA)
!     -----------------------------------
    theta(:) = 0.d0
!
    do i = 1, ndim
        do ino = 1, nnop
            theta(i) = theta(i)+ff(ino)*zr(ithet-1+ndim*(ino-1)+i)
        end do
!
!       ON REMPLACE E2 PAR THETA: REDRESSEMENT EN BORD DE FISSURE
        e2(i) = theta(i)
    end do
!
!    ON REORTHOGONALISE THETA
!    EN GARDANT SA NORME!
!
    call provec(e1, theta, temp)
    call normev(temp, norme)
    call normev(theta, norme_theta)
    call provec(temp, e1, theta)
    call normev(theta, norme)
    do i = 1, ndim
        theta(i) = norme_theta*theta(i)
    end do
!
!   NORMALISATION DE LA BASE
    call normev(e1, norme)
    call normev(e2, norme)
    call provec(e1, e2, e3)
    call normev(e3, norme)
!
!   OPTION 1 : ON REORTHOGONALISE EN CORRIGEANT N (DECOMMENTER)
!    call provec(e2, e3, e1)
!    call normev(e1, norme)
!
!   OPTION 2 : ON REORTHOGONALISE EN CORRIGEANT M (DECOMMENTER)
!    call provec(e3, e1, e2)
!    call normev(e2, norme)
!
!   OPTION 3 : DIRECTION TANGENTE = DIRECTION AU FOND
    do i = 1, ndim
        vec(i) = xg(i)-ptp(i)
    end do
!   Seule modification par rapport à avant :
!   on va prendre la projection dans le plan
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    cmp_hp = ddot(b_n, vec, b_incx, e3, b_incy)
    do i = 1, ndim
        vec(i) = vec(i)-cmp_hp*e3(i)
    end do
    call normev(vec, norme)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    sens = ddot(b_n, vec, b_incx, e2, b_incy)
    if (norme .ne. 0.d0) then
        if (sens .lt. 0.d0) then
            do i = 1, ndim
                e2(i) = -vec(i)
            end do
        else
            do i = 1, ndim
                e2(i) = vec(i)
            end do
        end if
    end if
!
!   ON ORTHOGONALISE
    call provec(e2, e3, e1)
    call normev(e1, norme)
!
! AJOUT CALCUL DE LA CONTRAINTE COHÉSIVE
!
    lambl(:) = 0.d0
    do ino = 1, nnop
        call indent(ino, ddls, ddlm, nnos, ii)
        do j = 1, ndim
            lambl(j) = lambl(j)+zr(idepl-1+ii+ddld+j)*ff(ino)
        end do
    end do
    call prmave(0, ptr, ndim, ndim, ndim,&
                lambl, ndim, lamb, ndim, ier)
!
!     ---------------------------------------------
!     3) CALCUL DU SAUT DE DEPLACEMENT
!     ---------------------------------------------
    grdep(:, :) = 0.d0
!
    do ino = 1, nnop
        am(:) = 0.d0
        call indent(ino, ddls, ddlm, nnos, ii)
        do j = 1, ndim
            do l = 1, ndim
                am(j) = am(j)+ptr(j, l)*zr(idepl-1+ddld+ndim+ii+l)
            end do
        end do
!
        do j = 1, ndim
            do l = 1, ndim
                grdep(j, l) = grdep(j, l)+dfdi(ino, l)*am(j)
            end do
        end do
    end do
!
!     -----------------------------------
!     4) CALCUL EFFECTIF DE G, K1, K2, K3
!     -----------------------------------
    g = 0.d0
    k1 = 0.d0
    k2 = 0.d0
    k3 = 0.d0
    g1 = 0.d0
    g2 = 0.d0
    g3 = 0.d0
    grde1 = 0.d0
    grde2 = 0.d0
    grde3 = 0.d0
    do j = 1, ndim
        do l = 1, ndim
            grde1 = grde1+e1(j)*grdep(j, l)*theta(l)
            grde2 = grde2+e2(j)*grdep(j, l)*theta(l)
            grde3 = grde3+e3(j)*grdep(j, l)*theta(l)
            g = g-lamb(j)*grdep(j, l)*theta(l)
        end do
    end do
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    lamb1 = ddot(b_n, lamb, b_incx, e1, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    lamb2 = ddot(b_n, lamb, b_incx, e2, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    lamb3 = ddot(b_n, lamb, b_incx, e3, b_incy)
    g1 = -lamb1*grde1
    g2 = -lamb2*grde2
    g3 = -lamb3*grde3
    gs2 = (lamb2*abs(grde2)-grde2*abs(lamb2))/2.d0
    gs3 = (lamb3*abs(grde3)-grde3*abs(lamb3))/2.d0
!
    if (ndim .eq. 3) then
        zr(igthet-1+1) = zr(igthet-1+1)+g*jac
        zr(igthet-1+2) = zr(igthet-1+2)+g1*jac
        zr(igthet-1+3) = zr(igthet-1+3)+gs2*jac
        zr(igthet-1+4) = zr(igthet-1+4)+gs3*jac
        zr(igthet-1+5) = zr(igthet-1+5)+g1*jac*coeff
        zr(igthet-1+6) = zr(igthet-1+6)+g2*jac*coeff
        zr(igthet-1+7) = zr(igthet-1+7)+g3*jac*coeff3
    else if (ndim .eq. 2) then
! PAS PROGRAMME POUR L INSTANT
!
        ASSERT(.false.)
!
    end if
end subroutine
