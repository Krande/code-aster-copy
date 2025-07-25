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

subroutine xside2(elrefp, ndim, coorse, elrese, igeom, &
                  he, nfh, ddlc, ddlm, nfe, &
                  basloc, nnop, npg, idecpg, typmod, &
                  imate, idepl, lsn, lst, &
                  nfiss, heavn, jstno, sig)
!
! aslint: disable=W1306,W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epstmc.h"
#include "asterfort/nbsigm.h"
#include "asterfort/rccoma.h"
#include "asterfort/reeref.h"
#include "asterfort/xcinem.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xkamat.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xnbddl.h"
    integer(kind=8) :: ndim, igeom, imate, nnop, npg, idepl, idecpg
    integer(kind=8) :: nfh, ddlc, nfe, nfiss, heavn(nnop, 5), jstno
    character(len=8) :: elrefp, elrese, typmod(*)
    real(kind=8) :: basloc(6*nnop), he(nfiss), coorse(*)
    real(kind=8) :: lsn(nnop), lst(nnop), sig(4, npg)
!
!
! person_in_charge: samuel.geniaut at edf.fr
!.......................................................................
!
!     BUT:  CALCUL DE L'OPTION SIEF_ELGA AVEC X-FEM EN 2D
!.......................................................................
! IN  ELREFP  : ELEMENT DE REFERENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNEES DES SOMMETS DU SOUS-ELEMENT
! IN  IGEOM   : COORDONNEES DES NOEUDS DE L'ELEMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ELT
! IN  NFH    : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  DDLM    : NOMBRE DE DDL PAR NOEUD MILIEU
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE AUX NOEUDS
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ELEMENT
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : MATERIAU CODE
! IN  IDEPL   : ADRESSE DU DEPLACEMENT A PARTIR DE LA CONF DE REF
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  NFISS   : NOMBRE DE FISSURES "VUES" PAR L'ÉLÉMENT
! IN  JHEAVN  : POINTEUR VERS LA DEFINITION HEAVISIDE
!
! OUT SIG     : CONTRAINTES (SIEF_ELGA)
!
!......................................................................
!
    character(len=2) :: k2bid
    character(len=16) :: phenom
    integer(kind=8) :: kpg, n, i, j, ino, iret, ipg, hea_se
    integer(kind=8) :: nno, nnos, npgbis, ddls, ddld, ddlm, ndimb
    integer(kind=8) :: jcoopg, jdfd2, jgano, idfde, ivf, ipoids, nbsig
    integer(kind=8) :: singu
    aster_logical :: axi
    real(kind=8) :: f(3, 3), eps(6)
    real(kind=8) :: instan, rac2
    real(kind=8) :: xg(ndim), xe(ndim), ff(nnop)
    real(kind=8) :: r8bi7(7), r8bi3(3)
    real(kind=8) :: dfdi(nnop, ndim)
    real(kind=8) :: fk(27, 3, 3), dkdgl(27, 3, 3, 3)
    real(kind=8) :: grad(3, 3)
    real(kind=8) :: zero, s, sth, d(4, 4), r, epsth(6)
    real(kind=8) :: ka, mu
    integer(kind=8) :: nnops
!
    data zero/0d0/
    data rac2/1.4142135623731d0/
!--------------------------------------------------------------------
!
!     ON AUTORISE UNIQUEMENT L'ISOTROPIE
    call rccoma(imate, 'ELAS', 1, phenom, iret)
    ASSERT(iret .eq. 0 .and. phenom .eq. 'ELAS')
!
!     INITIALISATIONS
    instan = r8vide()
    r8bi7(:) = 0.d0
    r8bi3(:) = 0.d0
!
!   NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, ddls, singu)
!
! ---- NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
!      -----------------------------------------
    nbsig = nbsigm()
!
    call elrefe_info(fami='RIGI', nnos=nnops)
!
! - INITIALISATION
    axi = typmod(1) .eq. 'AXIS'
!
!
    k2bid = ' '
!
    call elrefe_info(elrefe=elrese, fami='XINT', ndim=ndimb, nno=nno, nnos=nnos, &
                     npg=npgbis, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfde, &
                     jdfd2=jdfd2, jgano=jgano)
    ASSERT(npg .eq. npgbis .and. ndim .eq. ndimb)
!
! CALCUL DE L IDENTIFIANT DU SS ELEMENT
    hea_se = xcalc_code(nfiss, he_real=[he])
!
! - CALCUL POUR CHAQUE POINT DE GAUSS
!
    do kpg = 1, npg
!
!       NUMERO DU PG DANS LE FAMILLE 'XFEM'
        ipg = idecpg+kpg
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
!       CALCUL DES FF
        call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                    xe, ff, dfdi=dfdi)
!
!-----------------------------------------------------------------------
!         BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELT
!-----------------------------------------------------------------------
!
!
!         FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
        if (singu .gt. 0) then
            call xkamat(imate, ndim, axi, ka, mu)
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(1), &
                              lsn, lst, zr(igeom), ka, mu, ff, fk, dfdi, dkdgl)
        end if
!
!       CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE) ET DU DEPLACEMENT
!       RADIAL SI AXI (NECESSAIRE POUR LE CALCUL DES DEFORMATIONS EPS)
        r = 0.d0
        if (axi) then
            do ino = 1, nnop
                r = r+ff(ino)*zr(igeom-1+2*(ino-1)+1)
            end do
            ASSERT(r .ge. 0d0)
        end if
!
!       CALCUL DES DEFORMATIONS EPS
!
        call xcinem(axi, igeom, nnop, nnops, idepl, &
                    ndim, he, &
                    nfiss, nfh, singu, ddls, ddlm, &
                    fk, dkdgl, ff, dfdi, f, &
                    eps, grad, heavn)
!
!       CALCUL DES DEFORMATIONS THERMIQUES EPSTH
        epsth(:) = 0.d0
        call epstmc('XFEM', ndim, instan, '+', ipg, &
                    1, r8bi3, imate, 'CHAR_MECA_TEMP_R', &
                    epsth)
!
!       CALCUL DE LA MATRICE DE HOOKE (MATERIAU ISOTROPE)
        call dmatmc('XFEM', imate, instan, '+', ipg, &
                    1, r8bi3, nbsig, d)
!
!       VECTEUR DES CONTRAINTES
        do i = 1, nbsig
            s = zero
            sth = zero
            do j = 1, nbsig
                s = s+eps(j)*d(i, j)
                sth = sth+epsth(j)*d(i, j)
            end do
            sig(i, kpg) = s-sth
        end do
        sig(4, kpg) = sig(4, kpg)*rac2
!
    end do
!
end subroutine
