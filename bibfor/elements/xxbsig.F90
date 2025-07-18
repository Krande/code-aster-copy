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

subroutine xxbsig(elrefp, elrese, ndim, coorse, igeom, &
                  he, nfh, ddlc, ddlm, nfe, &
                  basloc, nnop, npg, sigma, &
                  lsn, lst, nfiss, heavn, jstno, &
                  codopt, ivectu, imate)
!
! aslint: disable=W1306,W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/indent.h"
#include "asterfort/lteatt.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xkamat.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/iimatu.h"
#include "asterfort/xnbddl.h"
    integer(kind=8) :: ndim, nfe, nfh, nfiss, nnop, npg
    integer(kind=8) :: ddlc, ddlm, heavn(nnop, 5), jstno
    integer(kind=8), optional :: imate
    integer(kind=8) :: codopt, igeom, ivectu
    real(kind=8) :: basloc(3*ndim*nnop), coorse(*), he(nfiss)
    real(kind=8) :: lsn(nnop), lst(nnop)
    real(kind=8) :: sigma(codopt*(2*ndim-1)+1, codopt*(npg-1)+1)
    character(len=8) :: elrefp, elrese
!
! person_in_charge: samuel.geniaut at edf.fr
!
!.......................................................................
!
!     BUT:  CALCUL  DU PRODUIT BT. SIGMA SUR UN SOUS-ELEMENT X-FEM
!.......................................................................
!
! IN  ELREFP  : ÉLÉMENT DE RÉFÉRENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNÉES DES SOMMETS DU SOUS-ÉLÉMENT
! IN  IGEOM   : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ÉLT
! IN  NFH     : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  DDLM    : NOMBRE DE DDL PAR NOEUD MILIEU (EN 2D)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ÉLÉMENT
! IN  SIGMA   : CONTRAINTES DE CAUCHY
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  CODOPT  : CODE DE L OPTION, POUR DIMENSIONNER LE TABLEAU SIGMA
!                 0 : REFE_FORC_NODA
!                 1 : FORC_NODA,CHAR_MECA_TEMP_R
!
! OUT IVECTU  : ADRESSE DU VECTEUR BT.SIGMA
!
!......................................................................
    integer(kind=8) :: kpg, i, ig, n, nn, m, dec(nnop)
    integer(kind=8) :: ddld, ddls, nno, nnops, nnos, npgbis, cpt
    integer(kind=8) :: idfde, ipoids, ivf, jcoopg, jdfd2, jgano, hea_se, i_dim
    integer(kind=8) :: singu, alp, ii
    real(kind=8) :: xg(ndim), xe(ndim), ff(nnop), jac
    real(kind=8) :: dfdi(nnop, ndim), f(3, 3)
    real(kind=8) :: def(6, nnop, ndim*(1+nfh+ndim)), voigt(2*ndim)
    real(kind=8) :: r
    real(kind=8) :: fk(27, 3, 3), dkdgl(27, 3, 3, 3), ka, mu
    aster_logical :: axi
!
    real(kind=8) :: rac2
    data rac2/1.4142135623731d0/
!--------------------------------------------------------------------
!
!     ATTENTION, EN 3D, ZR(VECTU) EST DIMENSIONNÉ DE
!     TELLE SORTE QU'ILS NE PRENNENT PAS EN COMPTE LES DDL SUR LES
!     NOEUDS MILIEU
!
!     NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, ddls, singu)
!
!     RECUPERATION DU NOMBRE DE NOEUDS SOMMETS DE L'ELEMENT PARENT
    call elrefe_info(fami='RIGI', nnos=nnops)
!
    if (ndim .eq. 2) then
        axi = lteatt('AXIS', 'OUI')
    else if (ndim .eq. 3) then
        axi = .false.
    end if
!     ADRESSE DES COORD DU SOUS ELT EN QUESTION
    call elrefe_info(elrefe=elrese, fami='XINT', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npgbis, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfde, &
                     jdfd2=jdfd2, jgano=jgano)
!
    ASSERT(npg .eq. npgbis)
    do n = 1, nnop
        call indent(n, ddls, ddlm, nnops, dec(n))
    end do
!
! CALCUL DE L IDENTIFIANT DU SS ELEMENT
    hea_se = xcalc_code(nfiss, he_real=[he])
!
!-----------------------------------------------------------------------
!     BOUCLE SUR LES POINTS DE GAUSS
    do kpg = 1, npg
!
!       COORDONNÉES DU PT DE GAUSS DANS LE REPÈRE RÉEL : XG
        xg(:) = 0.d0
        do i = 1, ndim
            do n = 1, nno
                xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+i)
            end do
        end do
!
!       JUSTE POUR CALCULER LES FF
!
        call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                    xe, ff, dfdi=dfdi)
!
!         FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
        if (nfe .gt. 0) then
            if (codopt .eq. 1) then
                ASSERT(present(imate))
                call xkamat(zi(imate), ndim, axi, ka, mu)
            else
                ka = 3.d0
                mu = sigma(1, 1)
            end if
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(1), &
                              lsn, lst, zr(igeom), ka, mu, ff, fk, dfdi, dkdgl)
        end if
!
! -     CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE)
        if (axi) then
            r = 0.d0
            do n = 1, nnop
                r = r+ff(n)*zr(igeom-1+2*(n-1)+1)
            end do
!
            ASSERT(r .gt. 0d0)
!          ATTENTION : LE POIDS N'EST PAS X R
!          CE SERA FAIT PLUS TARD AVEC JAC = JAC X R
        end if
!
!       COORDONNÉES DU POINT DE GAUSS DANS L'ÉLÉMENT DE RÉF PARENT : XE
!       ET CALCUL DE FF, DFDI, ET EPS
        f(:, :) = 0.d0
        do i = 1, 3
            f(i, i) = 1.d0
        end do
!
! - CALCUL DES ELEMENTS GEOMETRIQUES
!
!
!      CALCUL DES PRODUITS SYMETR. DE F PAR N,
        def(:, :, :) = 0.d0
        do n = 1, nnop
            cpt = 0
!         FONCTIONS DE FORME CLASSIQUES
            do i = 1, ndim
                cpt = cpt+1
                def(1, n, i) = f(i, 1)*dfdi(n, 1)
                def(2, n, i) = f(i, 2)*dfdi(n, 2)
                def(3, n, i) = 0.d0
                def(4, n, i) = (f(i, 1)*dfdi(n, 2)+f(i, 2)*dfdi(n, 1))/rac2
                if (ndim .eq. 3) then
                    def(3, n, i) = f(i, 3)*dfdi(n, 3)
                    def(5, n, i) = (f(i, 1)*dfdi(n, 3)+f(i, 3)*dfdi(n, 1))/rac2
                    def(6, n, i) = (f(i, 2)*dfdi(n, 3)+f(i, 3)*dfdi(n, 2))/rac2
                end if
            end do
!
!         TERME DE CORRECTION (3,3) AXI QUI PORTE EN FAIT SUR LE DDL 1
            if (axi) then
                def(3, n, 1) = f(3, 3)*ff(n)/r
            end if
!
!         ENRICHISSEMENT PAR HEAVYSIDE
            do ig = 1, nfh
                do i = 1, ndim
                    cpt = cpt+1
                    do m = 1, 2*ndim
                        def(m, n, cpt) = def(m, n, i)*xcalc_heav(heavn(n, ig), hea_se, heavn(n, 5))
                    end do
                    if (ndim .eq. 2) then
                        def(3, n, cpt) = 0.d0
                    end if
                end do
!
!   TERME DE CORRECTION (3,3) AXI PORTE SUR LE DDL 1+NDIM*IG
                if (axi) then
                    def(3, n, 1+ndim*ig) = f(3, 3)*ff(n)/r &
                                           *xcalc_heav(heavn(n, ig), hea_se, heavn(n, 5))
                end if
!
            end do
!
!         ENRICHISSEMENT PAR LES NFE FONTIONS SINGULIÈRES
            do alp = 1, ndim*nfe
                do i = 1, ndim
                    cpt = cpt+1
                    def(1, n, cpt) = f(i, 1)*dkdgl(n, alp, i, 1)
!
                    def(2, n, cpt) = f(i, 2)*dkdgl(n, alp, i, 2)
!
                    def(3, n, cpt) = 0.d0
!
                    def(4, n, cpt) = ( &
                                     f(i, 1)*dkdgl(n, alp, i, 2)+f(i, 2)*dkdgl(n, alp, i, 1) &
                                     )/rac2
!
                    if (ndim .eq. 3) then
                        def(3, n, cpt) = f(i, 3)*dkdgl(n, alp, i, 3)
                        def(5, n, cpt) = ( &
                                         f(i, 1)*dkdgl(n, alp, i, 3)+f(i, 3)*dkdgl(n, alp, i, 1) &
                                         )/rac2
                        def(6, n, cpt) = ( &
                                         f(i, 3)*dkdgl(n, alp, i, 2)+f(i, 2)*dkdgl(n, alp, i, 3) &
                                         )/rac2
                    end if
                end do
            end do
!
!   TERME DE CORRECTION (3,3) AXI PORTE SUR LE DDL 1+NDIM*(NFH+ALP)
!      EN AXI: ON PROJETTE L ENRICHISSEMENT VECTORIEL SUIVANT X
            if (axi) then
                do alp = 1, ndim*nfe
                    def(3, n, 1+ndim*(nfh+alp)) = f(3, 3)*fk(n, alp, 1)/r
                end do
            end if
!
            ASSERT(cpt .eq. ddld)
!
        end do
!
!       CALCULER LE JACOBIEN DE LA TRANSFO SSTET->SSTET REF
!       AVEC LES COORDONNEES DU SOUS-ELEMENT
        if (ndim .eq. 2) then
            call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                        jac)
        else if (ndim .eq. 3) then
            call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                        jac)
        end if
!
!       MODIFICATION DU JACOBIEN SI AXI
        if (axi) then
            jac = jac*r
        end if
!
        if (codopt .eq. 1) then
            do n = 1, 3
                voigt(n) = sigma(n, kpg)
            end do
            voigt(4) = sigma(4, kpg)*rac2
            if (ndim .eq. 3) then
                voigt(5) = sigma(5, kpg)*rac2
                voigt(6) = sigma(6, kpg)*rac2
            end if
        end if
!
        do n = 1, nnop
            nn = dec(n)
!
            do i = 1, ddld
                ii = iimatu(i, ndim, nfh, nfe)
                do m = 1, 2*ndim
                    if (codopt .eq. 1) then
                        zr(ivectu-1+nn+ii) = zr(ivectu-1+nn+ii)+def(m, &
                                                                    n, i)*voigt(m)*jac
                    else if (codopt .eq. 0) then
                        zr(ivectu-1+nn+ii) = zr(ivectu-1+nn+ii)+abs( &
                                             def(m, n, i)*sigma(1, 1)*jac)
                    else
                        ASSERT(.false.)
                    end if
                end do
!   POUR LES DDLS HEAVISIDE SEULEMENT ::
!    IL PEUT ARRIVER QUE L ESTIMATION DE LA FORCE DE REFERENCE TENDE VERS ZERO
!    ON REMPLACE LA VALEUR AU DDL HEAVISIDE I PAR LA VALEUR AU DDL PHYSIQUE I_DIM DU MEME NOEUD
                if (i .ge. (ndim+1) .and. i .le. (ndim+nfh*ndim) .and. codopt .eq. 0) then
                    i_dim = i-ndim*int((i-1)/ndim)
                    zr(ivectu-1+nn+i) = zr(ivectu-1+nn+i_dim)
                end if
            end do
!
        end do
!
!
    end do
!
end subroutine
