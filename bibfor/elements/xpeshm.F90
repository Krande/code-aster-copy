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
subroutine xpeshm(nno, nnop, nnops, ndim, nddls, &
                  nddlm, npg, igeom, jpintt, jpmilt, &
                  jheavn, ivf, ipoids, idfde, ivectu, &
                  ipesa, heavt, lonch, cnset, rho, &
                  axi, yaenrm, nfiss, nfh, jfisno)
!
! person_in_charge: daniele.colombo at ifpen.fr
!
! IN  ELREFP  : ÉLÉMENT DE RÉFÉRENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNÉES DES SOMMETS DU SOUS-ÉLÉMENT
! IN  IGEOM   : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ÉLT
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ÉLÉMENT
! IN  IVECTU  : INDICE DU SECONDE MEMBRE
! IN  NFISS   : NOMBRE DE FISSURES
! IN  NFH     : NOMBRE DE DDL HEAVISIDE
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elref1.h"
#include "asterfort/hmdeca.h"
#include "asterfort/indent.h"
#include "asterfort/reeref.h"
#include "asterfort/tecach.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
    aster_logical :: axi
    integer(kind=8) :: nse, ise, in, ino, nno, j, ndim, nfh, nfiss, jfisno
    integer(kind=8) :: nnop, nnops, n, nddls, nddlm, ipi, npg, ifiss
    integer(kind=8) :: igeom, jpintt, jpmilt, ivf, ipoids, idfde, jheavn
    integer(kind=8) :: ivectu, yaenrm, ifh, fisno(nnop, nfiss), iret, jtab(7)
    integer(kind=8) :: ipesa, dec1(nnop), dec2(nnop), icla, ienr, ncomp, dec
    integer(kind=8) :: heavt(*), lonch(10), cnset(*)
    integer(kind=8) :: heavn(nnop, 5), ig, ncompn, hea_se
    real(kind=8) :: xg(ndim), xe(ndim), coorse(30), dbid(nnop, ndim), he(nfiss)
    real(kind=8) :: ff(nnop), poids, rho, rx
    character(len=8) :: elrefp
!
    call elref1(elrefp)
!
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMEN
    nse = lonch(1)
!     RECUPERATION DE LA DEFINITION DES DDLS HEAVISIDES
    call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                itab=jtab)
    ncompn = jtab(2)/jtab(3)
    ASSERT(ncompn .eq. 5)
    do in = 1, nnop
        do ig = 1, ncompn
            heavn(in, ig) = zi(jheavn-1+ncompn*(in-1)+ig)
        end do
    end do
!
!     RECUPERATION DE LA CONNECTIVITÃ~I FISSURE - DDL HEAVISIDES
!     ATTENTION !!! FISNO PEUT ETRE SURDIMENTIONNÃ~I
    if (nfiss .eq. 1) then
        do ino = 1, nnop
            fisno(ino, 1) = 1
        end do
    else
        do ifh = 1, nfh
!    ON REMPLIT JUSQU'A NFH <= NFISS
            do ino = 1, nnop
                fisno(ino, ifh) = zi(jfisno-1+(ino-1)*nfh+ifh)
            end do
        end do
    end if
!
    call tecach('OOO', 'PHEAVTO', 'L', iret, nval=2, &
                itab=jtab)
    ncomp = jtab(2)
!
!     BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!     BOUCLE SUR LES 4/3 SOMMETS DU SOUS-TETRA/TRIA
        coorse(:) = 0.d0
        do in = 1, nno
            ino = cnset(nno*(ise-1)+in)
            do j = 1, ndim
                if (ino .lt. 1000) then
                    coorse(ndim*(in-1)+j) = zr(igeom-1+ndim*(ino-1)+j)
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    coorse(ndim*(in-1)+j) = zr(jpintt-1+ndim*(ino-1000- &
                                                              1)+j)
                else if (ino .gt. 2000 .and. ino .lt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-2000- &
                                                              1)+j)
                else if (ino .gt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-3000- &
                                                              1)+j)
                end if
            end do
        end do
!
!     DEFINITION DE LA FONCTION HEAVISIDE POUR CHAQUE SS-ELT
        do ifiss = 1, nfiss
            he(ifiss) = heavt(ncomp*(ifiss-1)+ise)
        end do
!
        hea_se = xcalc_code(nfiss, he_real=[he])
!
        do n = 1, nnop
            call indent(n, nddls, nddlm, nnops, dec1(n))
            call hmdeca(n, nddls, nddlm, nnops, dec2(n), &
                        dec)
        end do
! =====================================================================
! --- BOUCLE SUR LES POINTS D'INTEGRATION -----------------------------
! =====================================================================
        do ipi = 1, npg
!
!     COORDONNÉES DU PT DE GAUSS DANS LE REPÈRE RÉEL : XG
            xg(:) = 0.d0
            do j = 1, ndim
                do in = 1, nno
                    xg(j) = xg(j)+zr(ivf-1+nno*(ipi-1)+in)*coorse(ndim* &
                                                                  (in-1)+j)
                end do
            end do
!
!     XG -> XE (DANS LE REPERE DE l'ELREFP) ET VALEURS DES FF EN XE
            xe(:) = 0.d0
!
!     CALCUL DES FF ET DES DERIVEES DFDI POUR L'ELEMENT PARENTS
!     QUDRATIQUE (MECANIQUE)
            call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                        xe, ff, dbid)
!
!     CALCUL DU JACOBIEN DE LA TRANSFORMATION SS-ELT -> SS-ELT REF
            if (ndim .eq. 2) then
                call dfdm2d(nno, ipi, ipoids, idfde, coorse, &
                            poids)
            else if (ndim .eq. 3) then
                call dfdm3d(nno, ipi, ipoids, idfde, coorse, &
                            poids)
            end if
!
            poids = poids*rho*zr(ipesa)
            if (axi) then
                rx = 0.d0
                do ino = 1, nnop
                    rx = rx+zr(igeom+2*ino-2)*ff(ino)
                end do
                poids = poids*rx
!     TERMES CLASSIQUES
                do ino = 1, nnop
                    icla = dec1(ino)
                    ienr = dec2(ino)
!
                    zr(ivectu+icla+1) = zr(ivectu+icla+1)+poids*zr( &
                                        ipesa+2)*ff(ino)
!     TERMES HEAVISIDE
                    if (yaenrm .eq. 1) then
                        do ifh = 1, nfh
                            if (ino .le. nnops) then
                                zr(ivectu+ienr+1+(ndim+1)*ifh) = &
                                    zr(ivectu+ienr+1+(ndim+1)*ifh)+ &
                                    xcalc_heav(heavn(ino, ifh), hea_se, heavn(ino, 5))* &
                                    poids*zr(ipesa+2)*ff(ino)
                            else
                                zr(ivectu+ienr+1+ndim*ifh) = &
                                    zr(ivectu+ienr+1+ndim*ifh) &
                                    +xcalc_heav(heavn(ino, ifh), hea_se, heavn(ino, 5))* &
                                    poids*zr(ipesa+2)*ff(ino)
                            end if
                        end do
                    end if
                end do
            else
!     TERMES CLASSIQUES
                do ino = 1, nnop
                    icla = dec1(ino)
                    ienr = dec2(ino)
!
                    zr(ivectu+icla) = zr(ivectu+icla)+poids*zr(ipesa+1) &
                                      *ff(ino)
!
                    zr(ivectu+icla+1) = zr(ivectu+icla+1)+poids*zr( &
                                        ipesa+2)*ff(ino)
!
                    if (ndim .eq. 3) then
                        zr(ivectu+icla+2) = zr(ivectu+icla+2)+poids*zr( &
                                            ipesa+3)*ff(ino)
                    end if
!     TERMES HEAVISIDE
                    if (yaenrm .eq. 1) then
                        do ifh = 1, nfh
                            if (ino .le. nnops) then
                                zr(ivectu+ienr+(ndim+1)*ifh) = &
                                    zr(ivectu+ienr+(ndim+1)*ifh)+ &
                                    xcalc_heav(heavn(ino, ifh), hea_se, heavn(ino, 5))*poids* &
                                    zr(ipesa+1)*ff(ino)
                                zr(ivectu+ienr+1+(ndim+1)*ifh) = &
                                    zr(ivectu+ienr+1+(ndim+1)*ifh)+ &
                                    xcalc_heav(heavn(ino, ifh), hea_se, heavn(ino, 5))*poids* &
                                    zr(ipesa+2)*ff(ino)
                                if (ndim .eq. 3) then
                                    zr(ivectu+ienr+2+(ndim+1)*ifh) = &
                                        zr(ivectu+ienr+2+(ndim+1)*ifh)+ &
                                        xcalc_heav(heavn(ino, ifh), hea_se, heavn(ino, 5))*poids* &
                                        zr(ipesa+3)*ff(ino)
                                end if
                            else
                                zr(ivectu+ienr+ndim*ifh) = &
                                    zr(ivectu+ienr+ndim*ifh)+ &
                                    xcalc_heav(heavn(ino, ifh), hea_se, heavn(ino, 5))*poids* &
                                    zr(ipesa+1)*ff(ino)
!
                                zr(ivectu+ienr+ndim*ifh+1) = &
                                    zr(ivectu+ienr+ndim*ifh+1)+ &
                                    xcalc_heav(heavn(ino, ifh), hea_se, heavn(ino, 5))*poids* &
                                    zr(ipesa+2)*ff(ino)
                                if (ndim .eq. 3) then
                                    zr(ivectu+ienr+ndim*ifh+2) = &
                                        zr(ivectu+ienr+ndim*ifh+2)+ &
                                        xcalc_heav(heavn(ino, ifh), hea_se, heavn(ino, 5))*poids* &
                                        zr(ipesa+3)*ff(ino)
                                end if
                            end if
                        end do
                    end if
                end do
            end if
        end do
    end do
end subroutine
