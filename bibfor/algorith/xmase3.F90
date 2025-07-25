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

subroutine xmase3(elrefp, ndim, coorse, igeom, he, &
                  nfh, ddlc, nfe, basloc, nnop, &
                  npg, imate, lsn, lst, matuu, heavn, &
                  jstno, nnops, ddlm)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xkamat.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/iselli.h"
#include "asterfort/xnbddl.h"
#include "asterfort/iimatu.h"
#include "asterfort/indent.h"
    integer(kind=8) :: ndim, igeom, imate, nnop, npg, nfh, ddlc, nfe, heavn(27, 5)
    integer(kind=8) :: jstno, nnops, ddlm
    character(len=8) :: elrefp
    real(kind=8) :: basloc(9*nnop), he, coorse(*)
    real(kind=8) :: lsn(nnop), lst(nnop), matuu(*)
!
!     BUT:  CALCUL  DE L'OPTION MASS_MECA AVEC X-FEM EN 3D
!
!
! IN  ELREFP  : ÉLÉMENT DE RÉFÉRENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNÉES DES SOMMETS DU SOUS-ÉLÉMENT
! IN  IGEOM   : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ÉLT
! IN  DDLH    : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE AUX NOEUDS
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ÉLÉMENT
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : MATERIAU CODE
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  IDECPG  : POSITION DANS LA FAMILLE 'XFEM' DU 1ER POINT DE GAUSS
!               DU SOUS ELEMENT COURRANT (EN FAIT IDECPG+1)
! OUT MATUU   : MATRICE DE MASSE PROFIL
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: retour(1)
    integer(kind=8) :: kpg, kk, n, i, m, j, j1, kkd
    integer(kind=8) :: nno, nnos, npgbis, ddls, ddld, cpt, ndimb
    integer(kind=8) :: jcoopg, jdfd2, jgano, idfde, ivf, ipoids, hea_se
!
    real(kind=8) :: rho(1)
    real(kind=8) :: xg(ndim), xe(ndim), ff(nnop), jac
    real(kind=8) :: enr(ndim, nnop, 1+nfh+ndim*nfe)
    real(kind=8) :: fk(27, 3, 3), ka, mu
    integer(kind=8) :: alp, dec(nnop), nn, mn, ii, jj, irese, singu
    integer(kind=8) :: ddln, ij, kddl(ndim, 1+nfh+ndim*nfe)
!
    character(len=16) :: phenom
    character(len=8) :: elrese(6), fami(6)
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!--------------------------------------------------------------------
!
!     NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD SOMMET
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, ddls, singu)
    ddln = int(ddld/ndim)
    enr(:, :, :) = 0.d0
    kddl(:, :) = 0
!
! DECALAGES CALCULES EN AMONT: PERF
    do n = 1, nnop
        call indent(n, ddls, ddlm, nnops, dec(n))
    end do
!
! CALCUL DE L IDENTIFIANT DU SS ELEMENT
    hea_se = xcalc_code(1, he_real=[he])
!
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
!
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), &
                     ndim=ndimb, nno=nno, nnos=nnos, &
                     npg=npgbis, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfde, &
                     jdfd2=jdfd2, jgano=jgano)
!
    ASSERT(npg .eq. npgbis .and. ndim .eq. ndimb)
!
!
! - CALCUL POUR CHAQUE POINT DE GAUSS
    do kpg = 1, npg
!
!       COORDONNEES DU PT DE GAUSS DANS LE REPERE REEL : XG
        xg(:) = 0.d0
        do i = 1, ndim
            do n = 1, nno
                xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+ &
                                                             i)
            end do
        end do
!
!       JUSTE POUR CALCULER LES FF
        call reeref(elrefp, nnop, zr(igeom), xg, ndim, xe, ff)
!
        if (nfe .gt. 0) then
            call xkamat(imate, ndim, .false._1, ka, mu)
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he, &
                              lsn, lst, zr(igeom), ka, mu, ff, fk)
        end if
!
! - CALCUL DES ELEMENTS GEOMETRIQUES
!
!--------CALCUL DES FONCTIONS ENRICHIES--------------------------
        do n = 1, nnop
!         FONCTIONS DE FORME CLASSIQUES
            cpt = 0
            do i = 1, ndim
                cpt = cpt+1
                enr(i, n, 1) = ff(n)
                kddl(i, 1) = cpt
            end do
!         ENRICHISSEMENT PAR HEAVYSIDE
            do i = 1, ndim
                do j = 1, nfh
                    cpt = cpt+1
                    enr(i, n, 1+j) = ff(n)*xcalc_heav(heavn(n, j), hea_se, heavn(n, 5))
                    kddl(i, 1+j) = cpt
                end do
            end do
!         ENRICHISSEMENT PAR LES NFE FONTIONS SINGULIÈRES
            do alp = 1, nfe*ndim
                do i = 1, ndim
                    cpt = cpt+1
                    enr(i, n, 1+nfh+alp) = fk(n, alp, i)
                    kddl(i, 1+nfh+alp) = cpt
                end do
            end do
!
            ASSERT(cpt .eq. ddld)
!
        end do
!
!       POUR CALCULER LE JACOBIEN DE LA TRANSFO SSTET->SSTET REF
!       ON ENVOIE DFDM3D AVEC LES COORD DU SS-ELT
        call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                    jac)
!
!
!       ON RECUPERE LA MASSE VOLUMIQUE
!
        call rccoma(imate, 'ELAS', 1, phenom, retour(1))
        call rcvalb('RIGI', kpg, 1, '+', imate, &
                    ' ', phenom, 0, ' ', [0.d0], &
                    1, 'RHO', rho, retour, 1)
!
!
        do n = 1, nnop
            nn = dec(n)
            do ij = 1, ndim
                do i = 1, ddln
                    ii = iimatu(kddl(ij, i), ndim, nfh, nfe)
                    kkd = (nn+ii-1)*(nn+ii)/2
                    do j = 1, ddln
                        jj = iimatu(kddl(ij, j), ndim, nfh, nfe)
                        do m = 1, n
                            mn = dec(m)
                            if (m .eq. n) then
                                j1 = kddl(ij, i)
                            else
                                j1 = ddld
                            end if
                            if (jj .le. j1) then
                                kk = kkd+mn+jj
                                matuu(kk) = matuu(kk)+enr(ij, n, i)*enr(ij, m, j)*jac*rho(1)
                            end if
                        end do
                    end do
                end do
            end do
        end do
!
    end do
!
end subroutine
