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

subroutine nmasf3(nno, nbpg1, ipoids, ivf, idfde, &
                  imate, geom, deplm, sigm, vectu, &
                  compor)
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cast3d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elraga.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/invjac.h"
#include "asterfort/nmgeom.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nno, nbpg1, imate
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: ipoid2, ivf2, idfde2
    character(len=16) :: compor(*)
    real(kind=8) :: geom(3, nno)
    real(kind=8) :: deplm(3, nno), dfdi(nno, 3)
    real(kind=8) :: def(6, 3, nno)
    real(kind=8) :: sigm(78, nbpg1)
    real(kind=8) :: vectu(3, nno)
!.......................................................................
!
!     BUT:  CALCUL  DE L' OPTION FORC_NODA
!           EN HYPO-ELASTICITE EN 3D POUR LE HEXA8 SOUS INTEGRE
!           STABILITE PAR ASSUMED STRAIN
!.......................................................................
! IN  NNO     : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NBPG1   : NOMBRE DE POINTS DE GAUSS
! IN  POIDSG  : POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME
! IN  DFDE    : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  IMATE   : ADRESSE MATERIAU CODE
! IN  GEOM    : COORDONEES DES NOEUDS
! IN  DEPLM   : DEPLACEMENT A L'INSTANT PRECEDENT
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! OUT VECTU   : FORCES NODALES
!.......................................................................
!
!
    aster_logical :: grand, calbn, axi
    integer(kind=8) :: codre(1)
    character(len=16) :: nomres(2)
    character(len=32) :: phenom
    integer(kind=8) :: kpg, i, ii, ino, ia, j, k, kl, proj, nbpg2
    integer(kind=8) :: ndim, nnos, jgano, kp, iaa
    real(kind=8) :: d(6, 6), f(3, 3), eps(6), r, s
    real(kind=8) :: poids, poipg2(8)
    real(kind=8) :: jac, sigas(6, 8), invja(3, 3), bi(3, 8), hx(3, 4)
    real(kind=8) :: gam(4, 8), coopg2(24), h(8, 4), dh(4, 24)
    real(kind=8) :: qplus(72)
    real(kind=8) :: bn(6, 3, 8)
    real(kind=8) :: dfdx(8), dfdy(8), dfdz(8)
    real(kind=8) :: nu, nub, rac2, den
    real(kind=8) :: valres(2)
    data h/1.d0, 1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 1.d0, 1.d0,&
     &        1.d0, -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0,&
     &        1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0,&
     &       -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0, 1.d0, -1.d0/
!
! - INITIALISATION
!   ==============
!
!    PROJ : INDICATEUR DE LA PROJECTION
!           0 AUCUNE
!           1 ADS
!           2 ASQBI
!
    if (compor(1) .eq. 'ELAS            ') then
        proj = 2
    else
        proj = 1
    end if
    rac2 = sqrt(2.d0)
    grand = .false.
!
! - INITIALISATION HEXAS8
    call elraga('HE8', 'FPG8    ', ndim, nbpg2, coopg2, &
                poipg2)
    call elrefe_info(elrefe='HE8', fami='MASS', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=nbpg2, jpoids=ipoid2, jvf=ivf2, jdfde=idfde2, jgano=jgano)
!
! - CALCUL DES COEFFICIENTS BI (MOYENNE DES DERIVEES DES FCTS DE FORME)
    call r8inir(3*nno, 0.d0, bi, 1)
    den = 0.d0
    do kpg = 1, nbpg2
        call dfdm3d(nno, kpg, ipoid2, idfde2, geom, &
                    jac, dfdx, dfdy, dfdz)
        den = den+jac
        do ino = 1, nno
            bi(1, ino) = bi(1, ino)+jac*dfdx(ino)
            bi(2, ino) = bi(2, ino)+jac*dfdy(ino)
            bi(3, ino) = bi(3, ino)+jac*dfdz(ino)
        end do
    end do
    do i = 1, 3
        do ino = 1, nno
            bi(i, ino) = bi(i, ino)/den
        end do
    end do
!
! - CALCUL DES COEFFICIENTS GAMMA
!
    do i = 1, 4
        do k = 1, 3
            hx(k, i) = 0.d0
            do j = 1, nno
                hx(k, i) = hx(k, i)+h(j, i)*geom(k, j)
            end do
        end do
    end do
!
    do i = 1, 4
        do j = 1, nno
            s = 0.d0
            do k = 1, 3
                s = s+hx(k, i)*bi(k, j)
            end do
            gam(i, j) = 0.125d0*(h(j, i)-s)
        end do
    end do
!
! - CALCUL POUR LE POINT DE GAUSS CENTRAL
    kpg = 1
!
!  RECUP DU COEF DE POISSON POUR ASQBI
!
    call rccoma(imate, 'ELAS', 1, phenom, codre(1))
    nomres(1) = 'E'
    if (phenom .eq. 'ELAS') then
        nomres(2) = 'NU'
    else if (phenom .eq. 'ELAS_ISTR') then
        nomres(2) = 'NU_LT'
    else if (phenom .eq. 'ELAS_ORTH') then
        nomres(2) = 'NU_LT'
    else
        call utmess('F', 'ELEMENTS6_5', sk=phenom)
    end if
!
    call rcvalb('FPG1', 1, 1, '+', imate, &
                ' ', phenom, 0, ' ', [0.d0], &
                1, nomres(2), valres(2), codre, 1)
    if (codre(1) .eq. 0) then
        nu = valres(2)
    else
        call utmess('F', 'ELEMENTS4_72')
    end if
!
    nub = nu/(1.d0-nu)
!
    axi = .false.
    call nmgeom(3, nno, axi, grand, geom, &
                kpg, ipoids, ivf, idfde, deplm, &
                .true._1, poids, dfdi, f, eps, &
                r)
!
!      CALCUL DES PRODUITS SYMETR. DE F PAR N,
    do i = 1, nno
        do j = 1, 3
            def(1, j, i) = f(j, 1)*dfdi(i, 1)
            def(2, j, i) = f(j, 2)*dfdi(i, 2)
            def(3, j, i) = f(j, 3)*dfdi(i, 3)
            def(4, j, i) = (f(j, 1)*dfdi(i, 2)+f(j, 2)*dfdi(i, 1))/rac2
            def(5, j, i) = (f(j, 1)*dfdi(i, 3)+f(j, 3)*dfdi(i, 1))/rac2
            def(6, j, i) = (f(j, 2)*dfdi(i, 3)+f(j, 3)*dfdi(i, 2))/rac2
        end do
    end do
!
    do i = 1, 72
        qplus(i) = sigm(i+6, kpg)
    end do
!
    call r8inir(3*nno, 0.d0, vectu, 1)
    call r8inir(6*nbpg2, 0.d0, sigas, 1)
!
    calbn = .true.
!
!      OPERATEUR DE STABILISATION DU GRADIENT AUX 8 POINTS DE GAUSS
!
    do kpg = 1, nbpg2
        kp = 3*(kpg-1)
        call invjac(nno, kpg, ipoid2, idfde2, geom, &
                    invja, jac)
!
        do i = 1, 3
            dh(1, 3*(kpg-1)+i) = coopg2(3*kpg-1)*invja(3, i)+coopg2(3*kpg)*invja(2, i)
        end do
!
        do i = 1, 3
            dh(2, 3*(kpg-1)+i) = coopg2(3*kpg-2)*invja(3, i)+coopg2(3*kpg)*invja(1, i)
        end do
!
        do i = 1, 3
            dh(3, 3*(kpg-1)+i) = coopg2(3*kpg-2)*invja(2, i)+coopg2(3*kpg-1)*invja(1, i)
        end do
!
        do i = 1, 3
            dh(4, 3*(kpg-1)+i) = coopg2(3*kpg-2)*coopg2(3*kpg-1)*invja(3, i)+coopg2(3*kpg-1) &
                                &*coopg2(3*kpg)*invja(1, i)+coopg2(3*kpg-2)*coopg2(3*kpg)*&
                                & invja(2, i)
        end do
!
!
!  CALCUL DE BN AU POINT DE GAUSS KPG
!
        call cast3d(proj, gam, dh, def, nno, &
                    kpg, nub, nu, d, calbn, &
                    bn, jac, [0.d0])
!
!    CONTRAINTES DE HOURGLASS
!
        do i = 1, 6
            ii = 12*(i-1)
            do ia = 1, 4
                iaa = 3*(ia-1)
                do j = 1, 3
                    sigas(i, kpg) = sigas(i, kpg)+qplus(ii+iaa+j)*dh(ia, kp+j)
                end do
            end do
        end do
!
!     CALCUL DES FORCES INTERNES
!
        do i = 1, nno
            do j = 1, 3
                do kl = 1, 3
                    vectu(j, i) = vectu(j, i)+(def(kl, j, i)+bn(kl, j, i))*(sigas(kl, kpg)+sigm(kl,&
                                 &1))*jac+(rac2*def(kl+3, j, i)+bn(kl+3, j, i))*(sigas(kl+3, kpg)&
                                 &+sigm(kl+3, 1))*jac
                end do
            end do
        end do
!
    end do
!
end subroutine
