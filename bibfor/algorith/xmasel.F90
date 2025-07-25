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
subroutine xmasel(nnop, nfh, nfe, ddlc, igeom, &
                  imate, pintt, cnset, heavt, lonch, &
                  basloc, lsn, lst, matuu, heavn, &
                  jpmilt, jstno, nnops, ddlm)
    implicit none
#include "jeveux.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/xmase2.h"
#include "asterfort/xmase3.h"
    integer(kind=8) :: nnop, imate, igeom, jpmilt, jstno, ddlm, nnops
    integer(kind=8) :: nfh, nfe, ddlc, cnset(4*32), heavt(36), lonch(10), heavn(27, 5)
    real(kind=8) :: pintt(3*11), lsn(nnop)
    real(kind=8) :: lst(nnop), matuu(*), basloc(*)
!
!
!     BUT:  PRÉLIMINAIRES AU CALCUL DES OPTIONS RIGI_MECA_TANG,
!           RAPH_MECA ET FULL_MECA  EN HYPER-ELASTICITE AVEC X-FEM
!
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  IPOIDS  : POIDS DES POINTS DE GAUSS
! IN  IVF     : VALEUR  DES FONCTIONS DE FORME
! IN  nfh    : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  IGEOM   : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  NOMTE   : NOM DU TE
! IN  IMATE   : MATERIAU CODE
! IN  LGPG  : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!              CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  PINTT   : COORDONNÉES DES POINTS D'INTERSECTION
! IN  CNSET   : CONNECTIVITE DES SOUS-ELEMENTS
! IN  HEAVT   : VALEURS DE L'HEAVISIDE SUR LES SS-ELTS
! IN  LONCH   : LONGUEURS DES CHAMPS UTILISÉES
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! OUT MATUU   : MATRICE DE MASSE PROFIL
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: nse, npg, ndim
    integer(kind=8) :: j, ise, in, ino
    integer(kind=8) :: irese, nno
!
    real(kind=8) :: he, coorse(81)
!
    character(len=8) :: elrefp, elrese(6), fami(6)
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
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
!
!     ELEMENT DE REFERENCE PARENT : RECUP DE NDIM
    call elrefe_info(fami='RIGI', ndim=ndim)
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NPG
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), npg=npg, nno=nno)
!
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = lonch(1)
!
!       BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       BOUCLE SUR LES 4/3 SOMMETS DU SOUS-TETRA/TRIA
        do in = 1, nno
            ino = cnset((ndim+1)*(ise-1)+in)
            do j = 1, ndim
                if (ino .lt. 1000) then
                    coorse(ndim*(in-1)+j) = zr(igeom-1+ndim*(ino-1)+j)
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    coorse(ndim*(in-1)+j) = pintt(ndim*(ino-1000-1)+j)
                else if (ino .gt. 2000 .and. ino .lt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-2000-1)+j)
                else if (ino .gt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-3000-1)+j)
                end if
            end do
        end do
!
!       FONCTION HEAVYSIDE CSTE SUR LE SS-ELT
        he = heavt(ise)
!
!       DEBUT DE LA ZONE MEMOIRE DE SIG ET VI CORRESPONDANTE
!
        if (ndim .eq. 3) then
!
            call xmase3(elrefp, ndim, coorse, igeom, he, &
                        nfh, ddlc, nfe, basloc, nnop, &
                        npg, imate, lsn, lst, matuu, &
                        heavn, jstno, nnops, ddlm)
!
        else if (ndim .eq. 2) then
!
            call xmase2(elrefp, ndim, coorse, igeom, he, &
                        nfh, ddlc, nfe, basloc, nnop, &
                        npg, imate, lsn, lst, matuu, &
                        heavn, jstno, nnops, ddlm)
!
        end if
!
!
    end do
!
end subroutine
