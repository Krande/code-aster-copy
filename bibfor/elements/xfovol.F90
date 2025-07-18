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

subroutine xfovol(elrefp, ndim, coorse, igeom, he, &
                  ddlh, ddlc, singu, nnop, jlsn, &
                  jlst, heavn, iforc, itemps, ivectu, fonc, &
                  fono, imate, jbaslo, jstno)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/iselli.h"
#include "asterfort/lteatt.h"
#include "asterfort/reeref.h"
#include "asterfort/xkamat.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalc_code.h"
    character(len=8) :: elrefp
    real(kind=8) :: coorse(*)
    integer(kind=8) :: igeom, ndim, ddlh, ddlc, singu, nnop
    integer(kind=8) :: iforc, itemps, ivectu, jlsn, jlst, heavn(27, 5)
    integer(kind=8) :: imate, jbaslo, jstno
    real(kind=8) :: he
    aster_logical :: fonc, fono
!-----------------------------------------------------------------------
! FONCTION REALISEE : CALCUL DU SECOND MEMBRE AUX PG DU SOUS EL COURANT
!                     DANS LE CAS D'UNE FORCE VOLUMIQUE IMPOSEE SUR LES
!                     ELEMENTS X-FEM 2D ET 3D
!
!                    OPTIONS  :  CHAR_MECA_FR3D3D
!                                CHAR_MECA_FR2D2D
!                                CHAR_MECA_FF3D3D
!                                CHAR_MECA_FF2D2D
!
! IN  ELREFP  : ÉLÉMENT DE RÉFÉRENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNÉES DES SOMMETS DU SOUS-ÉLÉMENT
! IN  IGEOM   : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ÉLT
! IN  DDLH    : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  singu     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  JLSN    : INDICE DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  JLST    : INDICE DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN IFORC    : INDICE DE LA FORCE VOLUMIQUE
! IN ITEMPS   : INDICE DE L'INSTANT
! IN FONC     : .TRUE. SI LA FORCE EST FONCTION DE L'ESPACE OU DU TEMPS
! IN FONO     : .TRUE. SI LA FORCE EST FOURNIE AU NOEUD (.FALSE. AU PG)
! IN IVECTU   : INDICE DU SECONDE MEMBRE
!
!
!
    integer(kind=8) :: i, ino, ig, ier, j, n, mxstac
    integer(kind=8) :: ndimb, nno, nnos, nnops, npgbis, pos, irese, nfh
    integer(kind=8) :: jcoopg, ipoids, ivf, idfde, jdfd2, jgano, kpg, hea_se, alp
    real(kind=8) :: xe(ndim), xg(ndim), ff(nnop)
    real(kind=8) :: fk(27, 3, 3), ka, mu
    real(kind=8) :: forvol(ndim)
    real(kind=8) :: valpar(ndim+1), poids
    character(len=8) :: elrese(6), fami(6), nompar(ndim+1)
    aster_logical :: grdepl, axi
    parameter(mxstac=1000)
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!
!-----------------------------------------------------------------------
!     VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!     (VOIR CRS 1404)
    ASSERT(nnop .le. mxstac)
    ASSERT(ndim .le. mxstac)
!
    grdepl = .false.
!
    call elrefe_info(fami='RIGI', nnos=nnops)
!
    axi = lteatt('AXIS', 'OUI')
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
!     DEFINITION DE LA FONCTION HEAVISIDE POUR CHAQUE SS-ELT
    hea_se = xcalc_code(1, he_real=[he])
    nfh = ddlh/ndim
!
!     ------------------------------------------------------------------
!     BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELEMENT
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
!       CALCUL DES FONCTIONS D'ENRICHISSEMENT
!       -------------------------------------
!
        if (singu .gt. 0) then
            call xkamat(imate, ndim, axi, ka, mu)
            call xcalfev_wrap(ndim, nnop, zr(jbaslo), zi(jstno), he, &
                              zr(jlsn), zr(jlst), zr(igeom), ka, mu, ff, fk)
        end if
!
!       CALCUL DE LA FORCE VOLUMIQUE AU PG COURANT
!       ------------------------------------------
!
        forvol(:) = 0.d0
!
        if (fonc) then
!
!         FORCE AU PG COURANT A PARTIR DE LA FORCE FONCTION PAR ELEMENT
            do i = 1, ndim
                valpar(i) = xg(i)
            end do
            valpar(ndim+1) = zr(itemps)
            nompar(1) = 'X'
            nompar(2) = 'Y'
            if (ndim .eq. 3) then
                nompar(3) = 'Z'
                nompar(4) = 'INST'
                call fointe('FM', zk8(iforc), 4, nompar, valpar, &
                            forvol(1), ier)
                call fointe('FM', zk8(iforc+1), 4, nompar, valpar, &
                            forvol(2), ier)
                call fointe('FM', zk8(iforc+2), 4, nompar, valpar, &
                            forvol(3), ier)
            else if (ndim .eq. 2) then
                nompar(3) = 'INST'
                call fointe('FM', zk8(iforc), 3, nompar, valpar, &
                            forvol(1), ier)
                call fointe('FM', zk8(iforc+1), 3, nompar, valpar, &
                            forvol(2), ier)
            end if
!
        else
!
            if (fono) then
!           FORCE AU PG COURANT A PARTIR DE LA FORCE AUX NOEUDS
                do ino = 1, nnop
                    do j = 1, ndim
                        forvol(j) = forvol(j)+zr(iforc-1+ndim*(ino-1)+j) &
                                    *ff(ino)
                    end do
                end do
            else
!           FORCE FOURNIE AU PG
                do j = 1, ndim
                    forvol(j) = zr(iforc+j-1)
                end do
            end if
!
        end if
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
!
            end do
!
!         TERME HEAVISIDE
            do j = 1, ddlh
                pos = pos+1
                ig = j-nfh*int((j-1)/nfh)
                zr(ivectu-1+pos) = &
                    zr(ivectu-1+pos)+ &
                    xcalc_heav(heavn(ino, ig), hea_se, heavn(ino, 5))*forvol(j)*poids*ff(ino)
!
            end do
!
!         TERME SINGULIER
            do alp = 1, singu*ndim
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
