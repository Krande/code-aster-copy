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

subroutine xbsig(ndim, nnop, nfh, nfe, &
                 ddlc, ddlm, igeom, jpintt, &
                 cnset, heavt, lonch, basloc, sigma, &
                 nbsig, lsn, lst, ivectu, &
                 jpmilt, nfiss, jheavn, jstno, imate)
!
! aslint: disable=W1306,W1504
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/tecach.h"
#include "asterfort/xxbsig.h"
    integer(kind=8) :: ndim, nnop, nfh, nfe, ddlc, ddlm, igeom, nbsig, ivectu
    integer(kind=8) :: nfiss, jstno, imate
    integer(kind=8) :: cnset(4*32), heavt(*), lonch(10), jpintt, jpmilt, jheavn
    real(kind=8) :: basloc(*), sigma(*), lsn(nnop), lst(nnop)
!
!
! person_in_charge: samuel.geniaut at edf.fr
!
!      BSIGMC  -- CALCUL DES FORCES INTERNES B*SIGMA AUX NOEUDS
!                 DE L'ELEMENT DUES AU CHAMP DE CONTRAINTES SIGMA
!                 DEFINI AUX POINTS D'INTEGRATION DANS LE CADRE DE
!                 LA MÉTHODE X-FEM
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NFH     : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  DDLC    : NOMBRE DE DDLS DE CONTACT (PAR NOEUD)
! IN  DDLM    : NOMBRE DE DDL PAR NOEUD MILIEU
! IN  IGEOM   : COORDONEES DES NOEUDS
! IN  PINTT   : COORDONNÉES DES POINTS D'INTERSECTION
! IN  CNSET   : CONNECTIVITE DES SOUS-ELEMENTS
! IN  HEAVT   : VALEURS DE L'HEAVISIDE SUR LES SS-ELTS
! IN  LONCH   : LONGUEURS DES CHAMPS UTILISÉES
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE
! IN  SIGMA   : CONTRAINTES DE CAUCHY AUX POINTS DE GAUSS DES SOUS-ÉLTS
! IN  NBSIG   : NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  PMILT   : COORDONNEES DES POINTS MILIEUX
! IN  NFISS   : NOMBRE DE FISSURES "VUES" PAR L'ÉLÉMENT
!
! OUT IVECTU  : ADRESSE DU VECTEUR BT.SIGMA
!
!
!     VARIABLES LOCALES
    real(kind=8) :: he(nfiss), coorse(81)
    character(len=8) :: elrefp, elrese(6), fami(6)
    integer(kind=8) :: nse, idecpg, idebs, jtab(7), ncomp, iret
    integer(kind=8) :: ise, in, ino, npg, j, codopt
    integer(kind=8) :: ncompn, heavn(nnop, 5)
    integer(kind=8) :: irese, nno, ifiss, ig
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
!
    call elref1(elrefp)
!
!     NOMBRE DE COMPOSANTES DE PHEAVTO (DANS LE CATALOGUE)
    call tecach('OOO', 'PHEAVTO', 'L', iret, nval=2, &
                itab=jtab)
    ncomp = jtab(2)
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NNO ET NPG
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, npg=npg)
!
!   RECUPERATION DE LA DEFINITION DES FONCTIONS HEAVISIDES
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
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = lonch(1)
!
!       BOUCLE SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       BOUCLE SUR LES 4/3 SOMMETS DU SOUS-TETRA/TRIA
        do in = 1, nno
            ino = cnset(nno*(ise-1)+in)
!
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
!       FONCTION HEAVYSIDE CSTE POUR CHAQUE FISSURE SUR LE SS-ELT
        do ifiss = 1, nfiss
            he(ifiss) = heavt(ncomp*(ifiss-1)+ise)
        end do
!
!       DEBUT DE LA ZONE MÉMOIRE DE SIGMA CORRESPONDANTE
        idecpg = npg*(ise-1)
        idebs = nbsig*idecpg
        codopt = 1
        if (ndim .eq. 3) then
            ASSERT(nbsig .eq. 6)
        else if (ndim .eq. 2) then
            ASSERT(nbsig .eq. 4)
        end if
!
        call xxbsig(elrefp, elrese(ndim+irese), ndim, coorse, &
                    igeom, he, nfh, ddlc, ddlm, &
                    nfe, basloc, nnop, npg, sigma(idebs+1), &
                    lsn, lst, nfiss, &
                    heavn, jstno, codopt, ivectu, imate)
!
    end do
!
!
!.============================ FIN DE LA ROUTINE ======================
!
end subroutine
