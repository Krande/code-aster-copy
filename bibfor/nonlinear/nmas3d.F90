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
! aslint: disable=W1504
!
subroutine nmas3d(BEHInteg, &
                  fami, nno, npg, ipoids, ivf, &
                  idfde, geom, typmod, option, imate, &
                  compor, mult_comp, lgpg, carcri, instam, instap, &
                  deplm, deplp, angmas, sigm, vim, &
                  dfdi, def, sigp, vip, matuu, &
                  vectu, codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/caatdb.h"
#include "asterfort/calcdq.h"
#include "asterfort/cast3d.h"
#include "asterfort/codere.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elraga.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/invjac.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmgeom.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHinteg
    integer(kind=8) :: nno, imate, lgpg, codret, npg
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: ipoid2, ivf2, idfde2
    character(len=*) :: fami
    character(len=8) :: typmod(*)
    character(len=16) :: option
    character(len=16), intent(in) :: compor(*)
    character(len=16), intent(in) :: mult_comp
    real(kind=8), intent(in) :: carcri(*)
    real(kind=8) :: instam, instap
    real(kind=8) :: geom(3, nno)
    real(kind=8) :: deplm(3, nno), deplp(3, nno), dfdi(nno, 3)
    real(kind=8) :: def(6, 3, nno)
    real(kind=8) :: sigm(78, npg), sigp(78, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg)
    real(kind=8) :: matuu(*), vectu(3, nno), angmas(3)
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  CALCUL  DES OPTIONS RIGI_MECA_TANG, RAPH_MECA ET FULL_MECA
!           EN HYPO-ELASTICITE EN 3D POUR LE HEXA8 SOUS INTEGRE
!           STABILITE PAR ASSUMED STRAIN
!
! --------------------------------------------------------------------------------------------------
!
! IN  NNO     : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  POIDSG  : POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME
! IN  DFDE    : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  DFDK    : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  GEOM    : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  IMATE   : MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT PRECEDENT
! IN  INSTAP  : INSTANT DE CALCUL
! IN  DEPLM   : DEPLACEMENT A L'INSTANT PRECEDENT
! IN  DEPLP   : INCREMENT DE DEPLACEMENT
! IN  ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT DFDI    : DERIVEE DES FONCTIONS DE FORME  AU DERNIER PT DE GAUSS
! OUT DEF     : PRODUIT DER. FCT. FORME PAR F   AU DERNIER PT DE GAUSS
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
! OUT MATUU   : MATRICE DE RIGIDITE PROFIL (RIGI_MECA_TANG ET FULL_MECA)
! OUT VECTU   : FORCES NODALES (RAPH_MECA ET FULL_MECA)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    aster_logical :: grand, calbn, axi
    integer(kind=8) :: kpg, i, ii, ino, ia, j, k, kl, proj, cod(9), nbpg2
    integer(kind=8) :: nnos, jgano, kp, iaa, ndim2
    integer(kind=8), parameter :: ndimLdc = 3
    real(kind=8) :: d(6, 6), f(3, 3), eps(6), deps(6), r, s, sigma(6), sign(6)
    real(kind=8) :: poids, poipg2(8)
    real(kind=8) :: jac, sigas(6, 8), invja(3, 3), bi(3, 8), hx(3, 4)
    real(kind=8) :: gam(4, 8), coopg2(24), h(8, 4), dh(4, 24)
    real(kind=8) :: qplus(72), qmoins(72), dq(72)
    real(kind=8) :: bn(6, 3, 8)
    real(kind=8) :: pqx(4), pqy(4), pqz(4)
    real(kind=8) :: dfdx(8), dfdy(8), dfdz(8)
    real(kind=8) :: valres(2), nu, nub, den
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8) :: icodre(1)
    character(len=16) :: nomres(2)
    character(len=16) :: optios
    data h/1.d0, 1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 1.d0, 1.d0,&
     &        1.d0, -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0,&
     &        1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0,&
     &       -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0, 1.d0, -1.d0/
!
! --------------------------------------------------------------------------------------------------
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
    grand = .false.
    calbn = .false.

! - Prepare external state variables (geometry)
    call behaviourPrepESVAGeom(nno, npg, ndimLdc, &
                               ipoids, ivf, idfde, &
                               geom, BEHinteg, &
                               deplm, deplp)
!
! - INITIALISATION CODES RETOURS
    do kpg = 1, npg
        cod(kpg) = 0
    end do
!
! - INITIALISATION HEXAS8
    call elraga('HE8', 'FPG8    ', ndim2, nbpg2, coopg2, &
                poipg2)
    call elrefe_info(elrefe='HE8', fami='MASS', nno=nno, nnos=nnos, &
                     npg=nbpg2, jpoids=ipoid2, jvf=ivf2, jdfde=idfde2, jgano=jgano)
!
! - CALCUL DES COEFFICIENTS BI (MOYENNE DES DERIVEES DES FCTS DE FORME)
!
    bi(:, :) = 0.d0
!
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

! - CALCUL DES ELEMENTS GEOMETRIQUES
!     CALCUL DE DFDI,F,EPS,DEPS ET POIDS
!
    do j = 1, 6
        eps(j) = 0.d0
        deps(j) = 0.d0
    end do
    axi = .false.
    call nmgeom(3, nno, axi, grand, geom, &
                kpg, ipoids, ivf, idfde, deplm, &
                .true._1, poids, dfdi, f, eps, &
                r)
!
!     CALCUL DE DEPS
    call nmgeom(3, nno, axi, grand, geom, &
                kpg, ipoids, ivf, idfde, deplp, &
                .false._1, poids, dfdi, f, deps, &
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
    do i = 1, 3
        sign(i) = sigm(i, kpg)
    end do
    do i = 4, 6
        sign(i) = sigm(i, kpg)*rac2
    end do
!
! - LOI DE COMPORTEMENT
    if (option(1:9) .eq. 'RAPH_MECA') then
        optios = 'FULL_MECA'
    else
        optios = option
    end if

! - Set main parameters for behaviour (on point)
    call behaviourSetParaPoin(kpg, ksp, BEHinteg)

! - Integrator
    sigma = 0.d0
    call nmcomp(BEHinteg, &
                fami, kpg, ksp, ndimLdc, typmod, &
                imate, compor, carcri, instam, instap, &
                6, eps, deps, 6, sign, &
                vim(1, kpg), optios, angmas, &
                sigma, vip(1, kpg), 36, d, cod(kpg), mult_comp)
!
! - ERREUR D'INTEGRATION
    if (cod(kpg) .eq. 1) then
        goto 999
    end if
!
!  RECUP DU COEF DE POISSON POUR ASQBI
!
    if (proj .eq. 2) then
        nomres(1) = 'E'
        if (compor(1) .eq. 'ELAS') then
            nomres(2) = 'NU'
        else if (compor(1) .eq. 'ELAS_ISTR') then
            nomres(2) = 'NU_LT'
        else if (compor(1) .eq. 'ELAS_ORTH') then
            nomres(2) = 'NU_LT'
        else
            ASSERT(.false.)
        end if
!
!
        call rcvalb(fami, kpg, 1, '-', imate, &
                    ' ', compor(1), 0, ' ', [0.d0], &
                    1, nomres(2), valres(2), icodre, 1)
        if (icodre(1) .eq. 0) then
            nu = valres(2)
        else
            call utmess('F', 'ELEMENTS4_72')
        end if
!
        nub = nu/(1.d0-nu)
    end if
!
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
!
        call r8inir(300, 0.d0, matuu, 1)
!
!     CALCUL DE KC (MATRICE DE RIGIDITE AU CENTRE)
!     --------------------------------------------
        call caatdb(nno, def, d, def, poids, &
                    matuu)
!
!           CORRECTION DE LA MATRICE DE RIGIDITE
!                 CALCUL DE KSTAB
!     --------------------------------------------
!
!        CALCUL DES TERMES EVALUES AUX 8 POINTS DE GAUSS
        do kpg = 1, nbpg2
            call invjac(nno, kpg, ipoid2, idfde2, geom, &
                        invja, jac)
            do i = 1, 3
                dh(1, 3*(kpg-1)+i) = coopg2(3*kpg-1)*invja(3, i)+coopg2(3*kpg)*invja(2, i)
            end do
            do i = 1, 3
                dh(2, 3*(kpg-1)+i) = coopg2(3*kpg-2)*invja(3, i)+coopg2(3*kpg)*invja(1, i)
            end do
            do i = 1, 3
                dh(3, 3*(kpg-1)+i) = coopg2(3*kpg-2)*invja(2, i)+coopg2(3*kpg-1)*invja(1, i)
            end do
            do i = 1, 3
                dh(4, 3*(kpg-1)+i) = coopg2(3*kpg-2)*coopg2(3*kpg-1)*invja(3, i)+coopg2(3*kpg&
                                    &-1)*coopg2(3*kpg)*invja(1, i)+coopg2(3*kpg-2)*coopg2(&
                                    &3*kpg)*invja(2, i)
            end do
            call cast3d(proj, gam, dh, def, nno, &
                        kpg, nub, nu, d, calbn, &
                        bn, jac, matuu)
        end do
    end if
!
! - CALCUL DES FORCES INTERNES ET DES CONTRAINTES DE CAUCHY
!
    if (option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RAPH_MECA') then
!
!     INITIALISATION
        nbpg2 = 8
        do ia = 1, 4
            pqx(ia) = 0.d0
            pqy(ia) = 0.d0
            pqz(ia) = 0.d0
        end do
!
!     DEFORMATIONS GENERALISEES
        do ia = 1, 4
            do kl = 1, nno
                pqx(ia) = pqx(ia)+gam(ia, kl)*deplp(1, kl)
                pqy(ia) = pqy(ia)+gam(ia, kl)*deplp(2, kl)
                pqz(ia) = pqz(ia)+gam(ia, kl)*deplp(3, kl)
            end do
        end do
!
!      INCREMENT DES CONTRAINTES GENERALISEES
!
        call calcdq(proj, nub, nu, d, pqx, &
                    pqy, pqz, dq)
!
        do i = 1, 72
            qmoins(i) = sigm(i+6, 1)
            qplus(i) = qmoins(i)+dq(i)
        end do
!
        vectu(:, :) = 0.d0
        sigas(:, :) = 0.d0
!
        calbn = .true.
!
!      OPERATEUR DE STABILISATION DU GRADIENT AUX 8 POINTS DE GAUSS
!
        do kpg = 1, nbpg2
            kp = 3*(kpg-1)
            call invjac(nno, kpg, ipoid2, idfde2, geom, &
                        invja, jac)
            do i = 1, 3
                dh(1, 3*(kpg-1)+i) = coopg2(3*kpg-1)*invja(3, i)+coopg2(3*kpg)*invja(2, i)
            end do
            do i = 1, 3
                dh(2, 3*(kpg-1)+i) = coopg2(3*kpg-2)*invja(3, i)+coopg2(3*kpg)*invja(1, i)
            end do
            do i = 1, 3
                dh(3, 3*(kpg-1)+i) = coopg2(3*kpg-2)*invja(2, i)+coopg2(3*kpg-1)*invja(1, i)
            end do
            do i = 1, 3
                dh(4, 3*(kpg-1)+i) = coopg2(3*kpg-2)*coopg2(3*kpg-1)*invja(3, i)+coopg2(3*kpg&
                                    &-1)*coopg2(3*kpg)*invja(1, i)+coopg2(3*kpg-2)*coopg2(&
                                    &3*kpg)*invja(2, i)
            end do
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
                        vectu(j, i) = vectu(j, i)+(def(kl, j, i)+bn(kl, j, i))*(sigas(kl, kpg)+sigm&
                                     &a(kl))*jac+(rac2*def(kl+3, j, i)+bn(kl+3, j, i))*(sigas(kl&
                                     &+3, kpg)+sigma(kl+3)/rac2)*jac
                    end do
                end do
            end do
        end do
        do kl = 1, 3
            sigp(kl, 1) = sigma(kl)
            sigp(kl+3, 1) = sigma(kl+3)/rac2
        end do
        do i = 1, 72
            sigp(i+6, 1) = qplus(i)
        end do
    end if
!
999 continue
! - SYNTHESE DES CODES RETOURS
    call codere(cod, npg, codret)
!
end subroutine
