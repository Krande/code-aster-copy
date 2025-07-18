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
subroutine nmasf2(nno, npg, ipoids, ivf, idfde, &
                  geom, typmod, sigm, dfdi, vectu)
!
!
    implicit none
#include "asterf_types.h"
#include "asterfort/dfda2d.h"
#include "asterfort/iniqs4.h"
#include "asterfort/nmgeom.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: nno, npg
    character(len=8) :: typmod(*)
    integer(kind=8) :: ipoids, ivf, idfde
    real(kind=8) :: geom(2, nno)
    real(kind=8) :: dfdi(nno, 2)
    real(kind=8) :: sigm(10, npg)
    real(kind=8) :: vectu(2, nno)
!
!.......................................................................
!
!     BUT:  CALCUL  DES OPTIONS FORC_NODA
!           EN HYPO-ELASTICITE EN 2D POUR LE QUAD4 SOUS INTEGRE
!           STABILITE PAR ASSUMED STRAIN
!.......................................................................
! IN  NNO     : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IPOIDS  : INDICE POIDS DES POINTS DE GAUSS
! IN  IVFF     :INDICE VALEUR  DES FONCTIONS DE FORME
! IN  IDFDE  :INDICE DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  GEOM    : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODEELISATION
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! OUT DFDI    : DERIVEE DES FONCTIONS DE FORME  AU DERNIER PT DE GAUSS
! OUT VECTU   : FORCES NODALES
!.......................................................................
!
!
    aster_logical :: grand, axi
    integer(kind=8) :: kpg, n, i, j, kl, kpgs, proj, npgs
    real(kind=8) :: f(3, 3), eps(6), r
    real(kind=8) :: poids
    real(kind=8) :: rac2
!
!     AJ. VARIABLESPOIDSG
    real(kind=8) :: jac, sigas(4, 4), defc(4, 4, 2)
    real(kind=8) :: dh(8), gamma(8), coopg(8)
    real(kind=8) :: sdkdx(4), sdkdy(4), sdedx(4), sdedy(4), poi2sg(4)
    real(kind=8) :: sdfdy(4, 4), sdfdx(4, 4), sdfde(4, 4), sdfdk(4, 4)
    real(kind=8) :: qplus(6), defn(4, 4, 2), kron(3, 3), depbid(2, 4)
!
    data kron/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!
!
! - INITIALISATION
!   ==============
!
!    PROJ : INDICATEUR DE LA PROJECTION
!           0 AUCUNE
!           1 OPTIMAL BENDING
!           2 INCOMPRESSIBLE
    proj = 2
    rac2 = sqrt(2.d0)
    call r8inir(8, 0.d0, depbid, 1)
    grand = .false.
    axi = typmod(1) .eq. 'AXIS'
!
    do i = 1, 3
        do j = 1, 3
            f(i, j) = kron(i, j)
        end do
    end do
!
!
! - INITIALISATION QUAS4
    call iniqs4(nno, sdfde, sdfdk, poi2sg, coopg)
!
! - CALCUL DU VECTEUR GAMMA
    gamma(1) = ( &
               geom(1, 4)*(geom(2, 2)-geom(2, 3))+geom(1, 2)*(geom(2, 3)-geom(2, 4))+geom(1, 3)*(ge&
               &om(2, 4)-geom(2, 2)))/(2*(((geom(1, 4)-geom(1, 2))*(geom(2, 1)-geom(2, 3)))+(geom(&
               &1, 1)-geom(1, 3))*(geom(2, 2)-geom(2, 4))) &
               )
!
    gamma(2) = ( &
               geom(1, 4)*(geom(2, 3)-geom(2, 1))+geom(1, 3)*(geom(2, 1)-geom(2, 4))+geom(1, 1)*(ge&
               &om(2, 4)-geom(2, 3)))/(2*(((geom(1, 4)-geom(1, 2))*(geom(2, 1)-geom(2, 3)))+(geom(&
               &1, 1)-geom(1, 3))*(geom(2, 2)-geom(2, 4))) &
               )
!
    gamma(3) = ( &
               geom(1, 4)*(geom(2, 1)-geom(2, 2))+geom(1, 1)*(geom(2, 2)-geom(2, 4))+geom(1, 2)*(ge&
               &om(2, 4)-geom(2, 1)))/(2*(((geom(1, 4)-geom(1, 2))*(geom(2, 1)-geom(2, 3)))+(geom(&
               &1, 1)-geom(1, 3))*(geom(2, 2)-geom(2, 4))) &
               )
!
    gamma(4) = ( &
               geom(1, 3)*(geom(2, 1)-geom(2, 2))+geom(1, 1)*(geom(2, 2)-geom(2, 3))+geom(1, 2)*(ge&
               &om(2, 3)-geom(2, 1)))/(2*(((geom(1, 2)-geom(1, 4))*(geom(2, 1)-geom(2, 3)))-(geom(&
               &1, 1)-geom(1, 3))*(geom(2, 2)-geom(2, 4))) &
               )
!
!
! - CALCUL POUR LE POINT DE GAUSS CENTRAL
    kpg = 1
!
!
!
! - CALCUL DES ELEMENTS GEOMETRIQUES
!
!   CALCUL DE DFDI,R(EN AXI) ET POIDS
!
    call nmgeom(2, nno, axi, grand, geom, &
                kpg, ipoids, ivf, idfde, depbid, &
                .true._1, poids, dfdi, f, eps, &
                r)
!
!
!    OPERATEUR DE GRADIENT AU CENTRE
    do n = 1, nno
        do i = 1, 2
            defc(1, n, i) = f(i, 1)*dfdi(n, 1)
            defc(2, n, i) = f(i, 2)*dfdi(n, 2)
            defc(3, n, i) = 0.d0
            defc(4, n, i) = (f(i, 1)*dfdi(n, 2)+f(i, 2)*dfdi(n, 1))/rac2
        end do
    end do
!
!
!
! - CALCUL DE LA FORCE INTERIEURE ET DES CONTRAINTES DE CAUCHY
!
!    INITIALISATION
    npgs = 4
!
!    CONTRAINTES GENERALISEES
    do i = 1, 6
        qplus(i) = sigm(i+4, kpg)
    end do
!
!
!
!    OPERATEUR DE STABILISATION DU GRADIENT AU 4 POINTS DE GAUSS
    do kpgs = 1, npgs
!
!
        call dfda2d(kpgs, nno, poi2sg(kpgs), sdfde, sdfdk, &
                    sdedx, sdedy, sdkdx, sdkdy, sdfdx, &
                    sdfdy, geom, jac)
!
        dh(2*kpgs-1) = coopg(2*kpgs-1)*sdkdx(kpgs)+coopg(2*kpgs)*sdedx(kpgs)
        dh(2*kpgs) = coopg(2*kpgs-1)*sdkdy(kpgs)+coopg(2*kpgs)*sdedy(kpgs)
!
!
        do n = 1, nno
            do i = 1, 2
!
!         QUAS4 SANS PROJECTION
!         ---------------------
                if (proj .eq. 0) then
                    defn(1, n, i) = f(i, 1)*gamma(n)*dh(2*kpgs-1)
                    defn(2, n, i) = f(i, 2)*gamma(n)*dh(2*kpgs)
                    defn(3, n, i) = 0.d0
                    defn(4, n, i) = (f(i, 1)*gamma(n)*dh(2*kpgs)+f(i, 2)*gamma(n)*dh(2*kpgs-1))
!
!       OPTIMAL BENDING
!       ---------------
                else if (proj .eq. 1) then
                    defn(1, n, i) = f(i, 1)*gamma(n)*dh(2*kpgs-1)
                    defn(2, n, i) = f(i, 2)*gamma(n)*dh(2*kpgs)
                    defn(3, n, i) = 0.d0
                    defn(4, n, i) = 0.d0
!
!       INCOMPRESSIBLE
!       --------------
                else if (proj .eq. 2) then
                    defn(1, n, i) = f(i, 1)*gamma(n)*dh(2*kpgs-1)*(0.5d0)+f(i, 2)*gamma(n)*dh(2*k&
                                  &pgs)*(-0.5d0)
                    defn(2, n, i) = f(i, 2)*gamma(n)*dh(2*kpgs)*0.5d0+f(i, 1)*gamma(n)*dh(2*kpgs-1)&
                                   &*(-0.5d0)
                    defn(3, n, i) = 0.d0
                    defn(4, n, i) = 0.d0
                end if
!
            end do
        end do
!
!
!    CONTRAINTES DE HOURGLASS
!
!       QUAS4 SANS PROJECTION
!       ---------------------
        if (proj .eq. 0) then
            sigas(1, kpgs) = qplus(1)*dh(2*kpgs-1)+qplus(2)*dh(2*kpgs)
            sigas(2, kpgs) = qplus(3)*dh(2*kpgs-1)+qplus(4)*dh(2*kpgs)
            sigas(3, kpgs) = 0.d0
            sigas(4, kpgs) = (qplus(5)*dh(2*kpgs)+qplus(6)*dh(2*kpgs-1))/2
!
!       OPTIMAL BENDING
!       ---------------
        else if (proj .eq. 1) then
            sigas(1, kpgs) = qplus(1)*dh(2*kpgs-1)+qplus(2)*dh(2*kpgs)
            sigas(2, kpgs) = qplus(3)*dh(2*kpgs-1)+qplus(4)*dh(2*kpgs)
            sigas(3, kpgs) = 0.d0
            sigas(4, kpgs) = 0.d0
!
!       INCOMPRESSIBLE
!       --------------
        else if (proj .eq. 2) then
            sigas(1, kpgs) = (qplus(1)*dh(2*kpgs-1)+qplus(2)*dh(2*kpgs))
            sigas(2, kpgs) = (qplus(3)*dh(2*kpgs-1)+qplus(4)*dh(2*kpgs))
            sigas(3, kpgs) = 0.d0
            sigas(4, kpgs) = 0.d0
        end if
!
!
!
!   CALCUL DES FORCES INTERNES
!
        do n = 1, nno
            do i = 1, 2
                do kl = 1, 3
                    vectu(i, n) = vectu(i, n)+defc(kl, n, i)*sigas(kl, kpgs)*jac+defn(kl, n, i)*si&
                                 &gas(kl, kpgs)*jac
                end do
                vectu(i, n) = vectu(i, n)+defc(4, n, i)*sigas(4, kpgs)*jac*rac2+defn(4, n, i)*siga&
                             &s(4, kpgs)*jac
            end do
        end do
!
        do n = 1, nno
            do i = 1, 2
                do kl = 1, 3
                    vectu(i, n) = vectu(i, n)+defc(kl, n, i)*sigm(kl, kpg)*jac+defn(kl, n, i)*sigm(&
                                 &kl, kpg)*jac
                end do
                vectu(i, n) = vectu(i, n)+defc(4, n, i)*sigm(4, kpg)*rac2*jac+defn(4, n, i)*sigm(4,&
                             &kpg)*jac
            end do
        end do
    end do
!
!
end subroutine
