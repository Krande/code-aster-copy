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
subroutine xrigel(nnop, ddlh, nfe, ddlc, igeom, &
                  jpintt, cnset, heavt, lonch, basloc, &
                  lsn, lst, sig, matuu, jpmilt, &
                  heavn, jstno, imate)
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/nbsigm.h"
#include "asterfort/xrige2.h"
#include "asterfort/xrige3.h"
    integer(kind=8) :: nnop, igeom
    integer(kind=8) :: ddlh, nfe, ddlc, cnset(4*32), heavt(36), lonch(10)
    integer(kind=8) :: jpintt, jpmilt, heavn(27, 5), jstno, imate
    real(kind=8) :: lsn(nnop)
    real(kind=8) :: lst(nnop), matuu(*), sig(*), basloc(*)
!
!     BUT:  PRÉLIMINAIRES AU CALCUL DES OPTIONS RIGI_MECA_GE,
!           AVEC X-FEM
!
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  IPOIDS  : POIDS DES POINTS DE GAUSS
! IN  IVF     : VALEUR  DES FONCTIONS DE FORME
! IN  DDLH    : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  IGEOM   : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  NOMTE   : NOM DU TE
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  PINTT   : COORDONNÉES DES POINTS D'INTERSECTION
! IN  CNSET   : CONNECTIVITE DES SOUS-ELEMENTS
! IN  HEAVT   : VALEURS DE L'HEAVISIDE SUR LES SS-ELTS
! IN  LONCH   : LONGUEURS DES CHAMPS UTILISÉES
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  SIG     : CONTRAINTES DE CAUCHY
! OUT MATUU   : MATRICE DE MASSE PROFIL
!
!..............................................................
!
!
!
!
    integer(kind=8) :: nse, npg
    integer(kind=8) :: j, ise, in, ino, idebs
    integer(kind=8) :: idecpg, nbsig, ndim
    integer(kind=8) :: irese, nno
!
    real(kind=8) :: he, coorse(81)
!
    character(len=8) :: elrefp, elrese(6), fami(6)
!
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
! ----------------------------------------------------------------------
!
!     ATTENTION, DEPL ET VECTU SONT ICI DIMENSIONNÉS DE TELLE SORTE
!     QU'ILS NE PRENNENT PAS EN COMPTE LES DDL SUR LES NOEUDS MILIEU
!
    call elref1(elrefp)
!
!     ELEMENT DE REFERENCE PARENT : RECUP DE NDIM
    call elrefe_info(fami='RIGI', ndim=ndim)
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NPG
    if (.not. iselli(elrefp) .and. ndim .le. 2) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, npg=npg)
!
!     NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
    nbsig = nbsigm()
!
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = lonch(1)
!
!       BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       BOUCLE SUR LES 4/3 SOMMETS DU SOUS-TETRA/TRIA
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
!
!       FONCTION HEAVYSIDE CSTE SUR LE SS-ELT
        he = heavt(ise)
!
!       DEBUT DE LA ZONE MEMOIRE DE SIG CORRESPONDANTE
        idecpg = npg*(ise-1)
        idebs = nbsig*idecpg
!
        if (ndim .eq. 3) then
            ASSERT(nbsig .eq. 6)
            call xrige3(elrefp, ndim, coorse, igeom, he, &
                        heavn, ddlh, ddlc, nfe, basloc, &
                        nnop, npg, lsn, lst, sig(idebs+1), &
                        matuu, jstno, imate)
!
        else if (ndim .eq. 2) then
            ASSERT(nbsig .eq. 4)
            call xrige2(elrefp, elrese(ndim+irese), ndim, coorse, igeom, &
                        he, heavn, ddlh, ddlc, nfe, &
                        basloc, nnop, npg, lsn, lst, &
                        sig(idebs+1), matuu, jstno, imate)
!
        end if
!
    end do
!
end subroutine
