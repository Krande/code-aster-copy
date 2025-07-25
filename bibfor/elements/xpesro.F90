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

subroutine xpesro(elrefp, ndim, coorse, igeom, jheavt, ncomp, &
                  heavn, nfh, ddlc, nfe, nfiss, &
                  ise, nnop, jlsn, jlst, ivectu, &
                  fno, imate, jbaslo, jstno)
!
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/lteatt.h"
#include "asterfort/reeref.h"
#include "asterfort/xkamat.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
    character(len=8) :: elrefp
    real(kind=8) :: coorse(*)
    integer(kind=8) :: igeom, ndim, ddlc, nfe, nnop
    integer(kind=8) :: ivectu, jlsn, jlst, imate, jbaslo, jstno
    integer(kind=8) :: jheavt, nfh, nfiss, ise, heavn(27, 5), ncomp
    real(kind=8) :: fno(ndim*nnop)
!-----------------------------------------------------------------------
! FONCTION REALISEE : CALCUL DU SECOND MEMBRE AUX PG DU SOUS EL COURANT
!                     DANS LE CAS D'UNE FORCE VOLUMIQUE IMPOSEE SUR LES
!                     ELEMENTS X-FEM 2D ET 3D
!
!                          OPTIONS  :  CHAR_MECA_PESA_R
!                                      CHAR_MECA_ROTA_R
!
! IN  ELREFP  : ÉLÉMENT DE RÉFÉRENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNÉES DES SOMMETS DU SOUS-ÉLÉMENT
! IN  IGEOM   : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ÉLT
! IN  DDLH    : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ÉLÉMENT
! IN  JLSN    : INDICE DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  JLST    : INDICE DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  IDECPG  : POSITION DANS LA FAMILLE 'XFEM' DU 1ER POINT DE GAUSS
!               DU SOUS ELEMENT COURRANT (EN FAIT 1ER POINT : IDECPG+1)
! IN ITEMPS   : INDICE DE L'INSTANT
! IN FNO      : FORCES VOLUMIQUES AUX NOEUDS DE L'ELEMENT PARENT
! IN IVECTU   : INDICE DU SECONDE MEMBRE
!
!
    integer(kind=8) :: i, ino, ig, j, n
    integer(kind=8) :: ndimb, nno, nnos, nnops, npgbis, pos, ifiz, he(nfiss), hea_se
    integer(kind=8) :: jcoopg, ipoids, ivf, idfde, jdfd2, jgano, kpg
    real(kind=8) :: xe(ndim), xg(ndim), ff(nnop)
    real(kind=8) :: fk(27, 3, 3), ka, mu
    integer(kind=8) :: alp, singu
    real(kind=8) :: forvol(ndim)
    real(kind=8) :: poids, r
    character(len=8) :: elrese(6), fami(6)
    aster_logical :: grdepl, axi
    integer(kind=8) :: irese
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!-----------------------------------------------------------------------
    grdepl = .false.
!
    axi = lteatt('AXIS', 'OUI')
    singu = min(1, nfe)
!
    call elrefe_info(fami='RIGI', nnos=nnops)
!
!       FONCTION HEAVYSIDE CSTE POUR CHAQUE FISSURE SUR LE SS-ELT
    do ifiz = 1, nfiss
        he(ifiz) = zi(jheavt-1+ncomp*(ifiz-1)+ise)
    end do
    hea_se = xcalc_code(nfiss, he_inte=[he])
!
!     SOUS-ELEMENT DE REFERENCE
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), ndim=ndimb, nno=nno, &
                     nnos=nnos, npg=npgbis, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, &
                     jdfde=idfde, jdfd2=jdfd2, jgano=jgano)
    ASSERT(ndim .eq. ndimb)
!
!     ------------------------------------------------------------------
!     BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELEMENT COURANT
!     ------------------------------------------------------------------
!
    do kpg = 1, npgbis
!
!       COORDONNÉES DU PT DE GAUSS DANS LA CONFIG RÉELLE DU SE : XG
        xg(:) = 0.d0
        do i = 1, ndim
            do n = 1, nno
                xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+i)
            end do
        end do
!
!       CALCUL DES FF DE L'ELEMENT DE REFERENCE PARENT AU PG COURANT
        call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                    xe, ff)
!
!       POUR CALCULER LE JACOBIEN DE LA TRANSFO SS-ELT -> SS-ELT REF
!       ON ENVOIE DFDM3D/DFDM2D AVEC LES COORD DU SS-ELT
        if (ndim .eq. 3) call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                                     poids)
        if (ndim .eq. 2) call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                                     poids)
!
!
! -     CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE)
!       ET DU DEPL. RADIAL
        if (axi) then
            r = 0.d0
            do ino = 1, nnop
                r = r+ff(ino)*zr(igeom-1+2*(ino-1)+1)
            end do
!
!
!
            ASSERT(r .gt. 0d0)
        end if
!
        if (axi) then
            poids = poids*r
        end if
!
!       CALCUL DES FONCTIONS D'ENRICHISSEMENT
!       -------------------------------------
!
        if (nfe .gt. 0) then
            call xkamat(imate, ndim, axi, ka, mu)
            call xcalfev_wrap(ndim, nnop, zr(jbaslo), zi(jstno), real(he(1), 8), &
                              zr(jlsn), zr(jlst), zr(igeom), ka, mu, ff, fk)
        end if
!
!
!       CALCUL DE LA FORCE VOLUMIQUE AU PG COURANT
!       ------------------------------------------
!
        forvol(:) = 0.d0
        do ino = 1, nnop
            do j = 1, ndim
                forvol(j) = forvol(j)+fno(ndim*(ino-1)+j)*ff(ino)
            end do
        end do
!
!
!       CALCUL EFFECTIF DU SECOND MEMBRE
!       --------------------------------
!
        pos = 0
        do ino = 1, nnop
!
!         TERME CLASSIQUE
            do j = 1, ndim
                pos = pos+1
                zr(ivectu-1+pos) = zr(ivectu-1+pos)+forvol(j)*poids*ff(ino)
            end do
!
!         TERME HEAVISIDE
            do ig = 1, nfh
                do j = 1, ndim
                    pos = pos+1
                    zr(ivectu-1+pos) = zr(ivectu-1+pos)+ &
                                       xcalc_heav(heavn(ino, ig), hea_se, heavn(ino, 5)) &
                                       *forvol(j)*poids*ff(ino)
                end do
            end do
!
!         TERME SINGULIER
            do alp = 1, ndim*singu
                pos = pos+1
                do j = 1, ndim
                    zr(ivectu-1+pos) = zr(ivectu-1+pos)+fk(ino, alp, j)*forvol(j)*poids
                end do
            end do
!
!         ON SAUTE LES POSITIONS DES LAG DE CONTACT FROTTEMENT
            if (.not. iselli(elrefp)) then
                if (ino .le. nnops) pos = pos+ddlc
            else
                pos = pos+ddlc
            end if
!
!
        end do
!
!
    end do
!
!     ------------------------------------------------------------------
!     FIN BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELEMENT
!     ------------------------------------------------------------------
!
end subroutine
