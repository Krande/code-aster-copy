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
subroutine tanmat(alpha, beta, gamma, k1, k2, &
                  dmax1, dmax2, dam1, dam2, curv, &
                  dcurv, tanma2)
!
    implicit none
!
!     EVALUE L ENDOMMAGEMENT ET LA MATRICE TANGENTE POUR LE MODELE
!     DE COMPORTEMENT GLRC_DAMAGE
!
! IN  ALPHA : PARAMETRE  MATERIAU
! IN  BETA : PARAMETRE  MATERIAU
! IN  GAMMA : PARAMETRE  MATERIAU
! IN  K1 : SEUIL D ENDOMMAGEMENT SUP
! IN  K2 : SEUIL D ENDOMMAGEMENT INF
! IN  DMAX1 : ENDOMMAGEMNT MAXIMAL SUP
! IN  DMAX2 : ENDOMMAGEMNT MAXIMAL INF
! IN  CURV : TENSEUR DES COURBURES A L INSTANT T
! IN  DCURV : INCREMENT DU TENSEUR DES COURBURES
!
! IN/OUT DAM1 : ENDOMMAGEMENT SUP
! IN/OUT DAM2 : ENDOMMAGEMENT INF
!
! OUT TANMRP : MATRICE TANGENTE DANS LA BASE DES VECTEURS PROPRES
!
#include "asterfort/d2diag.h"
#include "asterfort/damage.h"
#include "asterfort/xifonc.h"
    integer(kind=8) :: i, j
!
    real(kind=8) :: tanma2(3, 3)
    real(kind=8) :: alpha, beta, gamma, k1, k2, dmax1, dmax2
    real(kind=8) :: dam1, dam2, curv(3), dcurv(3)
    real(kind=8) :: kappa(2, 2), kappvp(2), prp(2, 2), theta
    real(kind=8) :: tdm1rp(3, 3), tdm2rp(3, 3), trkapp, xvp1, xvp2
    real(kind=8) :: tanmrp(3, 3), aux, cc, ss, cs, prp33(3, 3)
    real(kind=8) :: kappn(2, 2), kappvn(2), thetn
    real(kind=8) :: tanmad(3, 3), tanmtn(3, 3), tanmto(3, 3)
    real(kind=8) :: cp(3, 3), cp2(3, 3), cp3(3, 3)
    real(kind=8) :: zerdam
!
    zerdam = 1.d-10
!
!     TENSEUR DES COURBURES A L INSTANT T : KAPPA ET T+1 KAPPN
!
    kappn(1, 1) = curv(1)+dcurv(1)
    kappn(2, 2) = curv(2)+dcurv(2)
    kappn(1, 2) = (curv(3)+dcurv(3))/2.d0
    kappn(2, 1) = kappn(1, 2)
!
    kappa(1, 1) = curv(1)
    kappa(2, 2) = curv(2)
    kappa(1, 2) = curv(3)/2.d0
    kappa(2, 1) = kappa(1, 2)
!
!     VALEUR PROPRE DES DEFORAMTIONS (PRP N EST PAS UTILISE)
    call d2diag(kappa, kappvp, prp, theta)
    call d2diag(kappn, kappvn, prp, thetn)
!
!     MISE A JOUR DE L ENDOMMAGEMENT DAM1 ET DAM2
!     CALCUL DES TERMES DE LA MATRICE TANGENTE VENANT DE L ENDO
!
    call damage(kappvp, +1, k1, dmax1, dam1, &
                tdm1rp, alpha, beta, gamma)
    call damage(kappvp, -1, k2, dmax2, dam2, &
                tdm2rp, alpha, beta, gamma)
!
!     MATRICE TANGENTE DANS LA BASE DES VECTEURS PROPRES
!
    trkapp = kappa(1, 1)+kappa(2, 2)
    xvp1 = xifonc(kappvp(1), dam1, dam2, gamma)
    xvp2 = xifonc(kappvp(2), dam1, dam2, gamma)
!
    tanmrp(1, 2) = 2*alpha*xifonc(trkapp, dam1, dam2, gamma)
    tanmrp(2, 1) = tanmrp(1, 2)
    tanmrp(1, 1) = tanmrp(1, 2)+2*beta*xvp1
    tanmrp(2, 2) = tanmrp(1, 2)+2*beta*xvp2
    tanmrp(1, 3) = 0.d0
    tanmrp(2, 3) = 0.d0
    tanmrp(3, 1) = 0.d0
    tanmrp(3, 2) = 0.d0
    aux = kappvp(2)-kappvp(1)
!
    if (abs(aux) .le. (zerdam*(abs(kappvp(1))+abs(kappvp(2))))) then
        tanmrp(3, 3) = 0.d0
    else
        tanmrp(3, 3) = beta*(xvp2*kappvp(2)-xvp1*kappvp(1))/aux
    end if
!
    cc = cos(theta)**2
    ss = sin(theta)**2
    cs = cos(theta)*sin(theta)
!
    prp33(1, 1) = cc
    prp33(2, 1) = ss
    prp33(3, 1) = -2.0d0*cs
    prp33(1, 2) = ss
    prp33(2, 2) = cc
    prp33(3, 2) = 2.0d0*cs
    prp33(1, 3) = cs
    prp33(2, 3) = -cs
    prp33(3, 3) = cc-ss
!
    if ((kappvp(1)*kappvn(1) .le. 0.d0) .or. (kappvp(2)*kappvn(2) .le. 0.d0)) then
        do j = 1, 3
            do i = 1, 3
                cp(i, j) = tdm1rp(i, j)+tdm2rp(i, j)
                cp2(i, j) = prp33(j, i)
            end do
        end do
!
        cp3 = matmul(cp, prp33)
        tanmad = matmul(cp2, cp3)
        cp3 = matmul(tanmrp, prp33)
        tanmto = matmul(cp2, cp3)
!
!
!     MATRICE TANGENTE DANS LA BASE DES VECTEURS PROPRES (NOUVELLE DEF)
!
        trkapp = kappn(1, 1)+kappn(2, 2)
        xvp1 = xifonc(kappvn(1), dam1, dam2, gamma)
        xvp2 = xifonc(kappvn(2), dam1, dam2, gamma)
!
        tanmrp(1, 2) = 2*alpha*xifonc(trkapp, dam1, dam2, gamma)
        tanmrp(2, 1) = tanmrp(1, 2)
        tanmrp(1, 1) = tanmrp(1, 2)+2*beta*xvp1
        tanmrp(2, 2) = tanmrp(1, 2)+2*beta*xvp2
        tanmrp(1, 3) = 0.d0
        tanmrp(2, 3) = 0.d0
        tanmrp(3, 1) = 0.d0
        tanmrp(3, 2) = 0.d0
        aux = kappvn(2)-kappvn(1)
!
        if (abs(aux) .le. zerdam*(abs(kappvn(1))+abs(kappvn(2)))) then
            tanmrp(3, 3) = 0.d0
        else
            tanmrp(3, 3) = beta*(xvp2*kappvn(2)-xvp1*kappvn(1))/aux
        end if
!
        cc = cos(thetn)**2
        ss = sin(thetn)**2
        cs = cos(thetn)*sin(thetn)
!
        prp33(1, 1) = cc
        prp33(2, 1) = ss
        prp33(3, 1) = -2.0d0*cs
        prp33(1, 2) = ss
        prp33(2, 2) = cc
        prp33(3, 2) = 2.0d0*cs
        prp33(1, 3) = cs
        prp33(2, 3) = -cs
        prp33(3, 3) = cc-ss
!
        do j = 1, 3
            do i = 1, 3
                cp2(i, j) = prp33(j, i)
            end do
        end do
!
        cp3 = matmul(tanmrp, prp33)
        tanmtn = matmul(cp2, cp3)
!
        do j = 1, 3
            do i = 1, 3
                tanma2(i, j) = (tanmto(i, j)+tanmtn(i, j))/2.d0+tanmad(i, j)
            end do
        end do
    else
        do j = 1, 3
            do i = 1, 3
                cp(i, j) = tanmrp(i, j)+tdm1rp(i, j)+tdm2rp(i, j)
                cp2(i, j) = prp33(j, i)
            end do
        end do
!
        cp3 = matmul(cp, prp33)
        tanma2 = matmul(cp2, cp3)
    end if
!
end subroutine
