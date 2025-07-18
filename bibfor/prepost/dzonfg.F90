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
subroutine dzonfg(nsommx, icnc, nelcom, numeli, inno, &
                  tbelzo, nbelzo, tbnozo, nbnozo)
!
!*********************************************
!              BUT DE CETTE ROUTINE :        *
! RECHERCHER LES EFS ET LES NOEUDS COMPOSANT *
! LES COUCHES 1, 2 ET 3                      *
!*********************************************
!
! IN  NSOMMX                 : NOMBRE DE SOMMETS MAX PAR EF
! IN  ICNC(NSOMMX+2,NELEM)   : CONNECTIVITE EF => NOEUDS CONNECTES
!     1ERE VALEUR = NBRE DE NOEUDS SOMMETS CONNECTES A L EF N°X
!     2EME VALEUR = 1 SI EF UTILE 0 SINON
!     SUITE EF N°X=>N° DE NOEUDS SOMMETS CONNECTES A X
!     EN 2D EF UTILE = QUAD OU TRIA
!     EN 3D EF UTILE = TETRA OU HEXA
! IN  NELCOM                 : NOMBRE MAX D'EF PAR NOEUD
! IN  NUMELI(NELCOM+2,NNOEM) : CONNECTIVITE NOEUD X=>EF CONNECTES A X
!     1ERE VALEUR = NBRE D EFS UTILES CONNECTES AU NOEUD N°X
!     2EME VALEUR = 0 NOEUD MILIEU OU NON CONNECTE A UN EF UTILE
!                   1 NOEUD SOMMET A L INTERIEUR + LIE A UN EF UTILE
!                   2 NOEUD SOMMET BORD + LIE A UN EF UTILE
!     CONNECTIVITE  NOEUD N°X=>N° DES EF UTILE CONNECTES A X
! IN  INNO                   : NOEUD CONSIDER
! OUT TBELZO(SOMME(NBELZO))  : TOUS LES EFS COMPOSANTS LA ZONE
! OUT NBELZO(3)              : NOMBRE D EFS DES COUCHES 1 2 ET 3
!                              NBELZO(1)=NBR D EFS COUCHE 1
!                              NBELZO(2)=NBR D EFS COUCHES 1 ET 2
!                              NBELZO(3)=NBR D EFS COUCHES 1, 2 ET 3
! OUT TBNOZO(SOMME(NBNOZO))  : TOUS LES NOEUDS COMPOSANTS LA ZONE
! OUT NBNOZO(3)              : NOMBRE DE NOEUDS DES COUCHES 1 2 ET 3
!
    implicit none
!
! DECLARATION GLOBALE
!
#include "asterf_types.h"
#include "asterfort/assert.h"
    integer(kind=8) :: nsommx, icnc(nsommx+2, *), nelcom, numeli(nelcom+2, *)
    integer(kind=8) :: inno
    integer(kind=8) :: tbelzo(1000), nbelzo(3), tbnozo(1000), nbnozo(3)
!
! DECLARATION LOCALE
!
    integer(kind=8) :: j, inel, nuef, ind, noeud, iel, elem, n
    integer(kind=8) :: nedep, nefin, nbelco, nbnoco
    aster_logical :: test
!
! 1 - EFS DES COUCHES 1, 2 ET 3
! 1.1 - COUCHE 1
!
    nbelzo(1) = numeli(1, inno)
    ASSERT(nbelzo(1) .le. 1000)
    do inel = 1, nbelzo(1)
        tbelzo(inel) = numeli(inel+2, inno)
    end do
!
! 1.2 - COUCHES 2 ET 3
!
    nedep = 1
    nefin = nbelzo(1)
!
    do j = 1, 2
        nbelco = 0
        do inel = nedep, nefin
            nuef = tbelzo(inel)
            do ind = 1, icnc(1, nuef)
                noeud = icnc(ind+2, nuef)
                do iel = 1, numeli(1, noeud)
                    elem = numeli(iel+2, noeud)
                    test = .true.
                    n = 1
60                  continue
                    if (n .le. (nefin+nbelco) .and. test) then
                        if (tbelzo(n) .eq. elem) test = .false.
                        n = n+1
                        goto 60
                    end if
                    if (test) then
                        nbelco = nbelco+1
                        ASSERT((nefin+nbelco) .le. 1000)
                        tbelzo(nefin+nbelco) = elem
                    end if
                end do
            end do
        end do
        if (j .eq. 1) then
            nbelzo(2) = nbelco
            nedep = nbelzo(1)+1
            nefin = nbelzo(1)+nbelzo(2)
        else
            nbelzo(3) = nbelco
        end if
    end do
!
! 2 - NOEUDS DES COUCHES 1 2 ET 3
!
    nbnozo(1) = 1
    nbnozo(2) = 0
    nbnozo(3) = 0
    tbnozo(1) = inno
    nbnoco = 1
!
    nedep = 0
!
    do j = 1, 3
        nefin = nedep+nbelzo(j)
        do inel = nedep+1, nefin
            nuef = tbelzo(inel)
            do ind = 1, icnc(1, nuef)
                noeud = icnc(ind+2, nuef)
                n = 1
                test = .true.
100             continue
                if (n .le. nbnoco .and. test) then
                    if (tbnozo(n) .eq. noeud) test = .false.
                    n = n+1
                    goto 100
                end if
                if (test) then
                    nbnozo(j) = nbnozo(j)+1
                    nbnoco = nbnoco+1
                    ASSERT(nbnoco .le. 1000)
                    tbnozo(nbnoco) = noeud
                end if
            end do
        end do
        nedep = nefin
    end do
!
    nbelzo(2) = nbelzo(1)+nbelzo(2)
    nbelzo(3) = nbelzo(2)+nbelzo(3)
end subroutine
