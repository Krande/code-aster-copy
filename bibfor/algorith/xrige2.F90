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

subroutine xrige2(elrefp, elrese, ndim, coorse, igeom, &
                  he, heavn, ddlh, ddlc, nfe, basloc, &
                  nnop, npg, lsn, lst, sig, &
                  matuu, jstno, imate)
!
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/lteatt.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xkamat.h"
#include "asterfort/xnbddl.h"
#include "asterfort/iimatu.h"
#include "asterfort/iipff.h"
    integer(kind=8) :: ndim, igeom, nnop, npg, ddlh, ddlc, nfe, heavn(27, 5)
    integer(kind=8) :: jstno, imate
    character(len=8) :: elrefp
    character(len=8) :: elrese
    real(kind=8) :: basloc(6*nnop), he, coorse(*)
    real(kind=8) :: lsn(nnop), lst(nnop), sig(48), matuu(*)
!
! person_in_charge: samuel.geniaut at edf.fr
!
!.......................................................................
!
!     BUT:  CALCUL  DE L'OPTION RIGI_GEOM AVEC X-FEM EN 2D
!.......................................................................
! IN  ELREFP  : ELEMENT DE REFERENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNEES DES SOMMETS DU SOUS-ELEMENT
! IN  IGEOM   : COORDONNEES DES NOEUDS DE L'ELEMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ELT
! IN  DDLH    : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE AUX NOEUDS
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ELEMENT
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  IDECPG  : POSITION DANS LA FAMILLE 'XFEM' DU 1ER POINT DE GAUSS
!               DU SOUS ELEMENT COURRANT (EN FAIT IDECPG+1)
! IN SIG     : CONTRAINTES DE CAUCHY
! OUT MATUU   : MATRICE DE MASSE PROFIL
!......................................................................
!
!
!
    integer(kind=8) :: kpg, kk, n, i, m, j, j1, kkd, ij
    integer(kind=8) :: nno, nnos, npgbis, ddls, ddld, ddldn, cpt, ndimb
    integer(kind=8) :: jcoopg, jdfd2, jgano, idfde, ivf, ipoids, hea_se, nfiss
    integer(kind=8) :: alp, ii, jj, k, nfh
    real(kind=8) :: tmp1
    real(kind=8) :: xg(ndim), xe(ndim), ff(nnop), jac
    real(kind=8) :: dfdi(nnop, ndim), pff(1+ddlh+nfe*ndim**2, nnop, ndim)
    real(kind=8) :: rac2
    real(kind=8) :: fk(27, 3, 3), dkdgl(27, 3, 3, 3), ka, mu
    integer(kind=8) :: nnops
!
    aster_logical :: axi
!
    data rac2/1.4142135623731d0/
!
!
!--------------------------------------------------------------------
!
!     NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD
    nfh = int(ddlh/ndim)
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, ddls, nfe)
    ddldn = 1+nfh+nfe*ndim**2
!     ELEMENT DE REFERENCE PARENT : RECUP DE NNOPS
    call elrefe_info(fami='RIGI', nnos=nnops)
!
    axi = lteatt('AXIS', 'OUI')
!
    call elrefe_info(elrefe=elrese, fami='XINT', ndim=ndimb, nno=nno, nnos=nnos, &
                     npg=npgbis, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfde, &
                     jdfd2=jdfd2, jgano=jgano)
    ASSERT(npg .eq. npgbis .and. ndim .eq. ndimb)
!
! CALCUL DE L IDENTIFIANT DU SS ELEMENT :: LE MULTI-HEAVISIDE N EST PAS PRIS EN COMPTE SANS FISNO
    nfiss = 1
    hea_se = xcalc_code(nfiss, he_real=[he])
!
!-----------------------------------------------------------------------
!     BOUCLE SUR LES POINTS DE GAUSS DU SOUS-TÉTRA
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
        call reeref(elrefp, nnop, zr(igeom), xg, ndim, xe, ff, dfdi=dfdi)
!
        if (nfe .gt. 0) then
!            call xkamat(imate, ndim, axi, ka, mu)
            imate = imate
            ka = 3.
            mu = sum(abs(sig(((kpg-1)*4+1):((kpg-1)*4+2))))/2.
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he, &
                              lsn, lst, zr(igeom), ka, mu, ff, fk, dfdi=dfdi, dkdgl=dkdgl)
        end if
!
!       POUR CALCULER LE JACOBIEN DE LA TRANSFO SSTET->SSTET REF
!       ON ENVOIE DFDM2D AVEC LES COORD DU SS-ELT
        call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                    jac)
!
!--------CALCUL DES FONCTIONS ENRICHIES--------------------------
!
        do i = 1, ndim
            do n = 1, nnop
!----------LA PARTIE CORRESPOINDANTE AU FF CLASSIQUES
                cpt = 1
                pff(cpt, n, i) = dfdi(n, i)
!----------LA PARTIE CORRESPONDANTE A L'ENRICHEMENT HEAVISIDE
                if (ddlh .eq. ndim) then
                    cpt = 2
                    pff(cpt, n, i) = dfdi(n, i)*xcalc_heav(heavn(n, 1), hea_se, heavn(n, 5))
                end if
!----------LA PARTIE CORRESPONDANTE A L'ENRICHEMENT CRACK-TIP
                do alp = 1, nfe*ndim
                    do k = 1, ndim
                        cpt = cpt+1
                        pff(cpt, n, i) = dkdgl(n, alp, i, k)
                    end do
                end do
                ASSERT(cpt .eq. ddldn)
            end do
        end do
!
        do n = 1, nnop
            do m = 1, n
                do i = 1, ddldn
                    do j = 1, ddldn
                        tmp1 = sig( &
                               (kpg-1)*4+1)*pff(i, n, 1)*pff(j, m, 1)+sig((kpg-1)*4+2)*pff(i, &
                                                                                   n, 2)*pff(j, m, &
                                                                                           2)+sig( &
                               (kpg-1)*4+4)*(pff(i, n, 1)*pff(j, m, 2)+pff(i, n, 2)*pff(j, m, 1) &
                                             )/rac2
!              STOCKAGE EN TENANT COMPTE DE LA SYMETRIE
                        if (m .eq. n) then
                            j1 = iipff(i, ndim, nfh, nfe)
                        else
                            j1 = ddldn
                        end if
                        if (iipff(j, ndim, nfh, nfe) .le. j1) then
                            do ij = 1, ndim
                                ii = iimatu((iipff(i, ndim, nfh, nfe)-1)*ndim+ij, ndim, nfh, nfe)
                                kkd = (ddls*(n-1)+ii-1)*(ddls*(n-1)+ii)/2
                                jj = iimatu((iipff(j, ndim, nfh, nfe)-1)*ndim+ij, ndim, nfh, nfe)
                                kk = kkd+ddls*(m-1)+jj
                                matuu(kk) = matuu(kk)+tmp1*jac
                            end do
                        end if
                    end do
                end do
            end do
        end do
!
    end do
!
end subroutine
