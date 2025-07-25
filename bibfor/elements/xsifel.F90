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
subroutine xsifel(elrefp, ndim, coorse, igeom, jheavt, &
                  ise, nfh, ddlc, ddlm, nfe, &
                  puls, basloc, nnop, idepl, lsn, &
                  lst, idecpg, igthet, fno, nfiss, &
                  jheavn, jstno)
!
! person_in_charge: samuel.geniaut at edf.fr
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/chauxi.h"
#include "asterfort/coor_cyl.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gbil3d.h"
#include "asterfort/gbilin.h"
#include "asterfort/indent.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/provec.h"
#include "asterfort/rcvad2.h"
#include "asterfort/rcvarc.h"
#include "asterfort/reeref.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xcinem.h"
#include "asterfort/xnbddl.h"
!
    character(len=8) :: elrefp
    integer(kind=8) :: igeom, ndim, nfh, ddlc, ddlm, nfe, nnop, idecpg, idepl, jheavn
    integer(kind=8) :: nfiss, jheavt, ise, jstno
    real(kind=8) :: fno(ndim*nnop), coorse(*)
    real(kind=8) :: basloc(3*ndim*nnop), lsn(nnop), lst(nnop), puls
!
!
!    - FONCTION REALISEE : CALCUL DU TAUX DE RESTITUTION D'ENERGIE
!                          ET DES FACTEURS D'INTENSITE DE CONTRAINTES
!                          PAR LA METHODE ENERGETIQUE G-THETA
!                          POUR LES ELEMENTS X-FEM (2D/3D)
!
! IN  ELREFP  : ELEMENT DE REFERENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNEEES DES SOMMETS DU SOUS-ELEMENT
! IN  IGEOM   : COORDONNEEES DES NOEUDS DE L'ELEMENT PARENT
! IN  NFH     : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  NFISS   : NOMBRE DE FISSURES "VUES" PAR L'ELEMENT
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIERES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE AUX NOEUDS
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  DEPL    : DEPLACEMENTS
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  IDECPG  : POSITION DANS LA FAMILLE 'XFEM' DU 1ER POINT DE GAUSS
!               DU SOUS ELEMENT COURRANT (EN FAIT 1ER POINT : IDECPG+1)
! IN  FNO     : FORCES VOLUMIQUES AUX NOEUDS DE L'ELEMENT PARENT
! OUT IGTHET  : G, K1, K2, K3
!
!
    integer(kind=8) :: ithet, imate, icomp, icour, igthet, jtab(7), ncomp
    integer(kind=8) :: ipoids, jcoopg, ivf, idfde, jdfd2, jgano, jsigse
    integer(kind=8) :: i, j, k, kpg, n, ino, iret, cpt, ig, ipg, in
    integer(kind=8) :: ndimb, nno, nnos, npgbis, ddld, ddls
    integer(kind=8) :: ifiss, isigi, ncmp
    integer(kind=8) :: heavn(nnop, 5), ncompn, hea_se
    integer(kind=8) :: iret1, iret2, iret3
    integer(kind=8) :: singu, alp
    real(kind=8) :: g, k1, k2, k3, coefk, coeff3, valres(4), alpha, he(nfiss)
    real(kind=8) :: devres(4), e, nu, lambda, mu, ka, c1, c2, c3, xg(ndim)
    real(kind=8) :: k3a
    real(kind=8) :: xe(ndim), ff(nnop), dfdi(nnop, ndim), f(3, 3)
    real(kind=8) :: eps(6), e1(3), e2(3), e3(3), p(3, 3)
    real(kind=8) :: invp(3, 3), rg, tg
    real(kind=8) :: courb(3, 3, 3), du1dm(3, 3)
    real(kind=8) :: du2dm(3, 3), sigin(6), dsigin(6, 3), epsref(6)
    real(kind=8) :: du3dm(3, 3), grad(ndim, ndim), dudm(3, 3), poids
    real(kind=8) :: dtdm(3, 3), tzero(3), dzero(3, 4), th
    real(kind=8) :: dudme(3, 4), dtdme(3, 4), du1dme(3, 4), du2dme(3, 4)
    real(kind=8) :: du3dme(3, 4), sigse(6*27), dfdx(27), dfdy(27), dfdz(27)
    real(kind=8) :: u1l(3), u2l(3), u3l(3), u1(3), u2(3), u3(3), r
    real(kind=8) :: depla(3), theta(3), tgudm(3), tpn(27), tref, tempg
    real(kind=8) :: ttrgu, ttrgv, dfdm(3, 4), cs, coef, rho, rac2
    real(kind=8) :: dtx, dty, dtz
    real(kind=8) :: fk(27, 3, 3), dkdgl(27, 3, 3, 3)
    integer(kind=8) :: icodre(4)
    character(len=16) :: nomres(4)
    character(len=8) :: elrese(6), fami(6), fami_se
    aster_logical :: lcour, axi, l_temp_noeu, l_not_zero
    integer(kind=8) :: irese, nnops, indenn, mxstac
    parameter(mxstac=1000)
!
    data nomres/'E', 'NU', 'ALPHA', 'RHO'/
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!     VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!     (VOIR CRS 1404)
!
    ASSERT(nnop .le. mxstac)
    ASSERT((3*ndim*nnop) .le. mxstac)
!
    rac2 = sqrt(2.d0)
!
!   ATTENTION, DEPL ET VECTU SONT ICI DIMENSIONNÉS DE TELLE SORTE
!   QU'ILS NE PRENNENT PAS EN COMPTE LES DDL SUR LES NOEUDS MILIEU
!
!   NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, &
                ddls, singu)
!
!   NOMBRE DE COMPOSANTES DE PHEAVTO (DANS LE CATALOGUE)
    call tecach('OOO', 'PHEAVTO', 'L', iret, nval=2, &
                itab=jtab)
    ncomp = jtab(2)
!
!   ELEMENT DE REFERENCE PARENT : RECUP DE NNOPS
    call elrefe_info(fami='RIGI', nnos=nnops)
!
    axi = lteatt('AXIS', 'OUI')
!
!   NOMBRE DE COMPOSANTES DES TENSEURS
    ncmp = 2*ndim
!
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
!
    call jevech('PTHETAR', 'L', ithet)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCOMPOR', 'L', icomp)
!
!   Verification du cadre theorique du calcul
    if (zk16(icomp-1+1) .ne. 'ELAS' .or. zk16(icomp-1+3) .ne. 'PETIT') then
        call utmess('F', 'RUPTURE1_24')
    end if
!
!   Sous-element de reference
    fami_se = fami(ndim+irese)
    if (nfe .gt. 0) then
        if (ndim .eq. 3 .and. count(zi((jstno-1+1):(jstno-1+nnop)) .eq. -2) .eq. 0) fami_se = &
            'XGEO'
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami_se, ndim=ndimb, nno=nno, nnos=nnos, &
                     npg=npgbis, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfde, &
                     jdfd2=jdfd2, jgano=jgano)
    ASSERT(ndim .eq. ndimb)
!
!   Recuperation de la contrainte initiale aux noeuds des sous-elts
    call tecach('ONO', 'PSIGISE', 'L', iret, iad=jsigse)
!
!   Indicateur de contrainte initiale
    isigi = 0
    if (jsigse .ne. 0) isigi = 1
!
    if (isigi .ne. 0) then
!       Passage de la contrainte initiale aux noeuds des sous-elts
!       dans un tableau local au sous-elt
        do i = 1, nno
            do j = 1, ncmp
                sigse(ncmp*(i-1)+j) = zr(jsigse-1+ncmp*nno*(ise-1)+ncmp*(i-1)+j)
            end do
        end do
!
    end if
!
!   TEMPERATURE DE REF
    call rcvarc(' ', 'TEMP', 'REF', 'XFEM', 1, &
                1, tref, iret)
    if (iret .ne. 0) tref = 0.d0
!
!   TEMPERATURE AUX NOEUDS PARENT
    l_temp_noeu = .false.
    do ino = 1, nnop
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, &
                    1, tpn(ino), iret)
        if (iret .ne. 0) tpn(ino) = 0.d0
    end do
    if (iret .eq. 0) l_temp_noeu = .true.
!
!   FONCTION HEAVYSIDE CSTE SUR LE SS-ÉLT ET PAR FISSURE
    do ifiss = 1, nfiss
        he(ifiss) = zi(jheavt-1+ncomp*(ifiss-1)+ise)
    end do
!
!   RECUPERATION DE LA DEFINITION DES FONCTIONS HEAVISIDE
    if (nfh .gt. 0) then
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
        ASSERT(ncompn .eq. 5)
        do ino = 1, nnop
            do ig = 1, ncompn
                heavn(ino, ig) = zi(jheavn-1+ncompn*(ino-1)+ig)
            end do
        end do
    end if
!
! CALCUL DE L IDENTIFIANT DU SS ELEMENT
    hea_se = xcalc_code(nfiss, he_real=[he])
!
!     ------------------------------------------------------------------
!     BOUCLE SUR LES POINTS DE GAUSS DU SOUS-TETRA
!     ------------------------------------------------------------------
!
    do kpg = 1, npgbis
!
!       INITIALISATIONS
        dtdm(:, :) = 0.d0
        du1dm(:, :) = 0.d0
        du2dm(:, :) = 0.d0
        du3dm(:, :) = 0.d0
        sigin(:) = 0.d0
        epsref(:) = 0.d0
        dsigin(:, :) = 0.d0
!
!
!       RECUPERATION DES DONNEES MATERIAUX
        ipg = idecpg+kpg
        call rcvad2('XFEM', ipg, 1, '+', zi(imate), &
                    'ELAS', 4, nomres, valres, devres, &
                    icodre)
        if (icodre(3) .ne. 0) then
            valres(3) = 0.d0
            devres(3) = 0.d0
        end if
        if (icodre(4) .ne. 0) then
!       on est sur que RHO est present en modal (te0297 : appel a cgverho)
            valres(4) = 0.d0
            devres(4) = 0.d0
        end if
        e = valres(1)
        nu = valres(2)
        alpha = valres(3)
        rho = valres(4)
        k3a = alpha*e/(1.d0-2.d0*nu)
!
        if (ndim .eq. 3 .or. (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) .or. axi) then
!
            lambda = nu*e/((1.d0+nu)*(1.d0-2.d0*nu))
            mu = e/(2.d0*(1.d0+nu))
            ka = 3.d0-4.d0*nu
            coefk = e/(1.d0-nu*nu)
            coeff3 = 2.d0*mu
            c1 = lambda+2.d0*mu
            c2 = lambda
            c3 = mu
            th = 1.d0
!
        else if (ndim .eq. 2 .and. lteatt('C_PLAN', 'OUI')) then
!
            ka = (3.d0-nu)/(1.d0+nu)
            mu = e/(2.d0*(1.d0+nu))
            coefk = e
            coeff3 = mu
            c1 = e/(1.d0-nu*nu)
            c2 = nu*c1
            c3 = mu
            th = (1.d0-2.d0*nu)/(1.d0-nu)
!
        end if
!
!       COORDONNEES DU PT DE GAUSS DANS LE REPERE REEL : XG
        xg(:) = 0.d0
        do i = 1, ndim
            do n = 1, nno
                xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+i)
            end do
        end do
!
!       CALCUL DES FF
        call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                    xe, ff, dfdi=dfdi)
!
!       POUR CALCULER LE JACOBIEN DE LA TRANSFO SS-ELT -> SS-ELT REF
!       AINSI QUE LES DERIVEES DES FONCTIONS DE FORMES DU SS-ELT
!       ON ENVOIE DFDM3D/DFDM2D AVEC LES COORD DU SS-ELT
        if (ndim .eq. 3) call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                                     poids, dfdx, dfdy, dfdz)
        if (ndim .eq. 2) call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                                     poids, dfdx, dfdy)
!
!       --------------------------------------
!       1) COORDONNÉES POLAIRES ET BASE LOCALE
!       --------------------------------------
!
!       FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
        if (singu .gt. 0) then
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(1), &
                              lsn, lst, zr(igeom), ka, mu, &
                              ff, fk, dfdi, dkdgl, elref=elrefp, &
                              kstop='C')
        end if
!
! -     CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE)
!       ET DU DEPL. RADIAL
        if (axi) then
            r = 0.d0
            do ino = 1, nnop
                r = r+ff(ino)*zr(igeom-1+2*(ino-1)+1)
            end do
!
            if (axi) then
                poids = poids*r
            end if
! Si R négative, on s'arrete
!
            ASSERT(r .gt. 0d0)
        end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    CALCUL DES COOR. CYL.
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        p(:, :) = 0.d0
        invp(:, :) = 0.d0
        call coor_cyl(ndim, nnop, basloc, zr(igeom), ff, &
                      p, invp, rg, tg, l_not_zero)
! BRICOLAGE POUR CALCULER LE SIGNE DE K2 QUAND NDIM=2
        if (ndim .eq. 2) then
            e1(:) = 0.d0
            e1(1:ndim) = p(1:ndim, 1)
            e2(:) = 0.d0
            e2(1:ndim) = p(1:ndim, 2)
            call provec(e1, e2, e3)
            p(3, 3) = e3(3)
            invp(3, 3) = e3(3)
        end if
!
!       ---------------------------------------------
!       2) CALCUL DU DEPLACEMENT ET DE SA DERIVEE (DUDM)
!       ---------------------------------------------
!
        depla(:) = 0.d0
!
!       CALCUL DE L'APPROXIMATION DU DEPLACEMENT
        do in = 1, nnop
!
            call indent(in, ddls, ddlm, nnops, indenn)
            cpt = 0
!           DDLS CLASSIQUES
            do i = 1, ndim
                cpt = cpt+1
                depla(i) = depla(i)+ff(in)*zr(idepl-1+indenn+cpt)
            end do
!           DDLS HEAVISIDE
            do ig = 1, nfh
                do i = 1, ndim
                    cpt = cpt+1
                    depla(i) = depla(i)+xcalc_heav(heavn(in, ig), hea_se, heavn(in, 5))*ff(in)*z&
                               &r(idepl-1+indenn+cpt)
                end do
            end do
!           DDL ENRICHIS EN FOND DE FISSURE
            do ig = 1, singu
                do alp = 1, ndim
                    cpt = cpt+1
                    do i = 1, ndim
                        depla(i) = depla(i)+fk(in, alp, i)*zr(idepl-1+indenn+cpt)
                    end do
                end do
            end do
!
        end do
!
!       CALCUL DU GRAD DE U AU POINT DE GAUSS
!
        call xcinem(axi, igeom, nnop, nnops, idepl, &
                    ndim, he, nfiss, nfh, singu, &
                    ddls, ddlm, fk, dkdgl, ff, &
                    dfdi, f, eps, grad, heavn)
!
!       ON RECOPIE GRAD DANS DUDM (CAR PB DE DIMENSIONNEMENT SI 2D)
        dudm(:, :) = 0.d0
        do i = 1, ndim
            do j = 1, ndim
                dudm(i, j) = grad(i, j)
            end do
        end do
!
!       ------------------------------------------------
!       3) CALCUL DU CHAMP THETA ET DE SA DERIVEE (DTDM)
!       ------------------------------------------------
!
        do i = 1, ndim
!
            theta(i) = 0.d0
            do ino = 1, nnop
                theta(i) = theta(i)+ff(ino)*zr(ithet-1+ndim*(ino-1)+i)
            end do
!
            do j = 1, ndim
                do ino = 1, nnop
                    dtdm(i, j) = dtdm(i, j)+zr(ithet-1+ndim*(ino-1)+i)*dfdi(ino, j)
                end do
            end do
        end do
!
!       --------------------------------------------------
!       4) CALCUL DU CHAMP DE TEMPERATURE ET DE SA DERIVEE
!       --------------------------------------------------
!
!       TEMPERATURE AU POINT DE GAUSS
        call rcvarc(' ', 'TEMP', '+', 'XFEM', ipg, &
                    1, tempg, iret)
        if (iret .ne. 0) tempg = 0.d0
        ttrgu = tempg-tref
        ttrgv = 0.d0
!
        tgudm(:) = 0.d0
!
!       cas de la varc TEMP, "continue" et donnee au noeud. Calcul
!       de ses derivees partielles
        if (l_temp_noeu) then
            do i = 1, ndim
                do ino = 1, nnop
                    tgudm(i) = tgudm(i)+dfdi(ino, i)*tpn(ino)
                end do
            end do
        end if
!
!       cas des varc DTX DTY DTZ, derivees partielles de la temperature
!       "discontinue". Ces varc sont donnees aux pg xfem
        call rcvarc(' ', 'DTX', '+', 'XFEM', ipg, &
                    1, dtx, iret1)
        if (iret1 .eq. 0) then
!           economisons les appels a rcvarc... si DTX est absent, pas
!           besoin de recuperer les autres composantes
            ASSERT(.not. l_temp_noeu)
            call rcvarc(' ', 'DTY', '+', 'XFEM', ipg, &
                        1, dty, iret2)
            ASSERT(iret2 .eq. 0)
            tgudm(1) = dtx
            tgudm(2) = dty
            if (ndim .eq. 3) then
                call rcvarc(' ', 'DTZ', '+', 'XFEM', ipg, &
                            1, dtz, iret3)
                ASSERT(iret3 .eq. 0)
                tgudm(3) = dtz
            end if
        end if
!
!       ------------------------------------------------
!       5) CALCUL DES CHAMPS AUXILIAIRES ET DE LEURS DERIVEES
!       -----------------------------------------------------
!
        if (ndim .eq. 2) then
!           NON PRISE EN COMPTE DE LA COURBURE
            lcour = .false.
        else if (ndim .eq. 3) then
!           PRISE EN COMPTE DE LA COURBURE
            lcour = .true.
!           RECUPERATION DU TENSEUR DE COURBURE
            call jevech('PCOURB', 'L', icour)
            do i = 1, ndim
                do j = 1, ndim
                    courb(i, 1, j) = zr(icour-1+ndim*(i-1)+j)
                    courb(i, 2, j) = zr(icour-1+ndim*(i+3-1)+j)
                    courb(i, 3, j) = zr(icour-1+ndim*(i+6-1)+j)
                end do
            end do
        end if
!
        call chauxi(ndim, mu, ka, rg, tg, &
                    invp, lcour, courb, du1dm, du2dm, &
                    du3dm, u1l, u2l, u3l)
!
!       CHAMPS SINGULIERS DANS LA BASE GLOBALE
        u1(:) = 0.d0
        u2(:) = 0.d0
        u3(:) = 0.d0
        do i = 1, ndim
            do j = 1, ndim
                u1(i) = u1(i)+p(i, j)*u1l(j)
                u2(i) = u2(i)+p(i, j)*u2l(j)
                if (ndim .eq. 3) u3(i) = u3(i)+p(i, j)*u3l(j)
            end do
        end do
!
!       --------------------------------------------------
!       6) CALCUL DES TERMES LIEES A LA CONTRAINTE INITIALE
!       --------------------------------------------------
!
        if (isigi .ne. 0) then
!
!           Calcul de la contrainte initiale (somme sur les noeuds du ss-elt)
            do i = 1, nno
                do j = 1, ncmp
                    sigin(j) = sigin(j)+sigse(ncmp*(i-1)+j)*zr(ivf-1+nno*(kpg-1)+i)
                end do
            end do
!
!           Calcul du gradient de sigma initial (somme sur les noeuds du ss-elt)
            do i = 1, nno
                do j = 1, ncmp
                    dsigin(j, 1) = dsigin(j, 1)+sigse(ncmp*(i-1)+j)*dfdx(i)
                    dsigin(j, 2) = dsigin(j, 2)+sigse(ncmp*(i-1)+j)*dfdy(i)
                    if (ndim .eq. 3) dsigin(j, 3) = dsigin(j, 3)+sigse(ncmp*(i-1)+j)*dfdz(i)
                end do
            end do
!
!           Traitements particuliers des termes croises
            do i = 4, ncmp
                sigin(i) = sigin(i)*rac2
                do j = 1, ndim
                    dsigin(i, j) = dsigin(i, j)*rac2
                end do
            end do
!
            epsref(1) = -(1.d0/e)*(sigin(1)-(nu*(sigin(2)+sigin(3))))
            epsref(2) = -(1.d0/e)*(sigin(2)-(nu*(sigin(3)+sigin(1))))
            epsref(3) = -(1.d0/e)*(sigin(3)-(nu*(sigin(1)+sigin(2))))
            epsref(4) = -(1.d0/mu)*sigin(4)
            if (ndim .eq. 3) then
                epsref(5) = -(1.d0/mu)*sigin(5)
                epsref(6) = -(1.d0/mu)*sigin(6)
            end if
!
        end if
!
!       -----------------------------------------------------------
!       7) CALCUL DES FORCES VOLUMIQUES ET DE LEURS DERIVEES (DFDM)
!       -----------------------------------------------------------
!
        dfdm(:, :) = 0.d0
        do ino = 1, nnop
            do j = 1, ndim
                do k = 1, ndim
                    dfdm(j, k) = dfdm(j, k)+fno(ndim*(ino-1)+j)*dfdi(ino, k)
                end do
!               VALEUR DE LA FORCE DANS LA QUATRIEME COLONNE :
                dfdm(j, 4) = dfdm(j, 4)+fno(ndim*(ino-1)+j)*ff(ino)
            end do
        end do
!
        if (axi) then
            dfdm(3, 3) = dfdm(1, 4)/r
        end if
!
!       ---------------------------------------------
!       8) CALCUL DE G, K1, K2, K3 AU POINT DE GAUSS
!       --------------------------------------------
!
        tzero(:) = 0.d0
        dzero(:, :) = 0.d0
!
!       POUR L'APPEL A GIL3D/GBIL2D, ON STOCKE LES CHAMPS
!       EN DERNIERE COLONNE DES MATRICES DES DERIVEES DES CHAMPS
!       ON ETEND LES MATRICES :
!       EX : DUDM DE DIM (3,3) -> DUDME DE DIM (3,4)
        dudme(:, :) = 0.d0
        dtdme(:, :) = 0.d0
        du1dme(:, :) = 0.d0
        du2dme(:, :) = 0.d0
        du3dme(:, :) = 0.d0
        do i = 1, ndim
            do j = 1, ndim
                dudme(i, j) = dudm(i, j)
                dtdme(i, j) = dtdm(i, j)
                du1dme(i, j) = du1dm(i, j)
                du2dme(i, j) = du2dm(i, j)
                du3dme(i, j) = du3dm(i, j)
            end do
            dudme(i, 4) = depla(i)
            dtdme(i, 4) = theta(i)
            du1dme(i, 4) = u1(i)
            du2dme(i, 4) = u2(i)
            du3dme(i, 4) = u3(i)
        end do
!
        if (axi) then
            dudme(3, 3) = dudme(1, 4)/r
            dtdme(3, 3) = dtdme(1, 4)/r
            du1dme(3, 3) = du1dme(1, 4)/r
            du2dme(3, 3) = du2dme(1, 4)/r
            du3dme(3, 3) = du3dme(1, 4)/r
        end if
!
        if (ndim .eq. 3) then
!
            coef = 2.d0
            call gbil3d(dudme, dudme, dtdme, dfdm, dfdm, &
                        tgudm, tgudm, ttrgu, ttrgu, poids, &
                        sigin, dsigin, epsref, c1, c2, &
                        c3, k3a, alpha, coef, rho, &
                        puls, g)
            zr(igthet) = zr(igthet)+g
!
            coef = 1.d0
            call gbil3d(dudme, du1dme, dtdme, dfdm, dzero, &
                        tgudm, tzero, ttrgu, ttrgv, poids, &
                        sigin, dsigin, epsref, c1, c2, &
                        c3, k3a, alpha, coef, rho, &
                        puls, k1)
            zr(igthet+4) = zr(igthet+4)+k1*coefk
            zr(igthet+1) = zr(igthet+1)+k1*sqrt(coefk)
!
            coef = 1.d0
            call gbil3d(dudme, du2dme, dtdme, dfdm, dzero, &
                        tgudm, tzero, ttrgu, ttrgv, poids, &
                        sigin, dsigin, epsref, c1, c2, &
                        c3, k3a, alpha, coef, rho, &
                        puls, k2)
            zr(igthet+5) = zr(igthet+5)+k2*coefk
            zr(igthet+2) = zr(igthet+2)+k2*sqrt(coefk)
!
            coef = 1.d0
            call gbil3d(dudme, du3dme, dtdme, dfdm, dzero, &
                        tgudm, tzero, ttrgu, ttrgv, poids, &
                        sigin, dsigin, epsref, c1, c2, &
                        c3, k3a, alpha, coef, rho, &
                        puls, k3)
            zr(igthet+6) = zr(igthet+6)+k3*coeff3
            zr(igthet+3) = zr(igthet+3)+k3*sqrt(coeff3)
!
        else if (ndim .eq. 2) then
!
!           POUR G, COEF = 2
            coef = 2.d0
            cs = 1.d0
            call gbilin('XFEM', ipg, zi(imate), dudme, dudme, &
                        dtdme, dfdm, tgudm, poids, sigin, &
                        dsigin, epsref, c1, c2, c3, &
                        cs, th, coef, rho, puls, &
                        axi, g)
!
!           POUR K1, COEF = 1
            coef = 1.d0
            cs = 0.5d0
            call gbilin('XFEM', ipg, zi(imate), dudme, du1dme, &
                        dtdme, dfdm, tgudm, poids, sigin, &
                        dsigin, epsref, c1, c2, c3, &
                        cs, th, coef, rho, puls, &
                        axi, k1)
            k1 = k1*coefk
!
!           POUR K2, COEF = 1
            coef = 1.d0
            cs = 0.5d0
            call gbilin('XFEM', ipg, zi(imate), dudme, du2dme, &
                        dtdme, dfdm, tgudm, poids, sigin, &
                        dsigin, epsref, c1, c2, c3, &
                        cs, th, coef, rho, puls, &
                        axi, k2)
            k2 = k2*coefk
            if (e3(3) .lt. 0) k2 = -k2
!
            zr(igthet) = zr(igthet)+g
            zr(igthet+1) = zr(igthet+1)+k1/sqrt(coefk)
            zr(igthet+2) = zr(igthet+2)+k2/sqrt(coefk)
            zr(igthet+3) = zr(igthet+3)+k1
            zr(igthet+4) = zr(igthet+4)+k2
!
        end if
!
    end do
!
!     ------------------------------------------------------------------
!     FIN DE LA BOUCLE SUR LES POINTS DE GAUSS DU SOUS-TETRA
!     ------------------------------------------------------------------
!
end subroutine
