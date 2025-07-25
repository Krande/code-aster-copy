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
subroutine dfort2(nsommx, icnc, noeu1, tbelzo, nbelt, &
                  tbnozo, nbnozo, nbnoe, xy, aire, &
                  energi, pe)
!
!********************************************************************
!                  BUT DE CETTE ROUTINE :                           *
!       CALCULER LE DEGRE DE LA SINGULARITE PE AU NOEUD NOEU1       *
! 1) CALCUL DE L ENERGIE MOYENNE ENER SUR UN CERCLE DE CENTRE NOEU1 *
!    ET POUR DIFFERENTS RAYONS RAYZ                                 *
! 2) COMPARAISON DE CETTE ENERGIE CALCULEE AVEC L ENERGIE THEORIQUE *
!    EN POINTE DE FISSURE ENER_THEO=K*(RAYZ**(2*(PE-1)))+C          *
!    ET CALCUL DU COEFFICENT PE PAR LA METHODE DES MOINDRES CARRES  *
!********************************************************************
!
! IN  NSOMMX               : NOMBRE DE SOMMETS MAX PAR EF
! IN  ICNC(NSOMMX+2,NELEM) : CONNECTIVITE EF => NOEUDS CONNECTES
!     1ERE VALEUR = NBRE DE NOEUDS SOMMETS CONNECTES A L EF N°X
!     2EME VALEUR = 1 SI EF UTILE 0 SINON
!     CONNECTIVITE  EF N
!     EN 2D EF UTILE = QUAD OU TRIA
!     EN 3D EF UTILE = TETRA OU HEXA
! IN  NDIM                 : DIMENSION DU PROBLEME
! IN  NOEU1                : NOEUD CONSIDERE
! IN  TBELZO(NBELT)        : EFS DES COUCHES 1, 2 ET 3
! IN  NBELT                : NBR D EFS DES COUCHES 1,2,3
! IN  TBNOZO(NBNOE)        : NOEUDS DES COUCHES 1, 2 ET 3
! IN  NBNOZO(3)            : NBR DE NOEUDS DES COUCHES 1, 2 ET 3
! IN  XY(3,NNOEM)          : COORDONNEES DES NOEUDS
! IN  AIRE(NELEM)          : SURFACE DE CHAQUE EF
! IN  ENERGI(NELEM)        : ENERGIE SUR CHAQUE EF
! OUT PE                   : DEGRE DE LA SINGULARITE
!
    implicit none
!
! DECLARATION GLOBALE
!
#include "asterfort/assert.h"
#include "asterfort/dcalph.h"
#include "asterfort/dcqpri.h"
#include "asterfort/dcspri.h"
#include "asterfort/dinter.h"
    integer(kind=8) :: nsommx, icnc(nsommx+2, *), noeu1
    integer(kind=8) :: nbnozo(3), nbelt, nbnoe
    integer(kind=8) :: tbelzo(nbelt), tbnozo(nbnoe)
    real(kind=8) :: xy(3, *), aire(*), energi(*), pe
!
! DECLARATION LOCALE
!
    integer(kind=8) :: i, j, inno, iint, inel, nuef, noeu2, nedep, nefin
    integer(kind=8) :: nsomm, nbint
    integer(kind=8) :: ipoi1, ipoi2, ipoi4, nint, ip1
    parameter(nbint=10)
    real(kind=8) :: coord(2), coor(2, 4), coorin(2, 2)
    real(kind=8) :: delta(3), dist, rayz(nbint), ray
    real(kind=8) :: airtot, sprim
    real(kind=8) :: ener(nbint), epsir
!---------------------------------------------------------------------------------
    epsir = 1.d0+1.d-15
!
! 1 - COORDONNEES DU NOEUD CONSIDERE INNO
!
    do i = 1, 2
        coord(i) = xy(i, noeu1)
    end do
!
! 2 - CALCUL DES RAYONS DES COUCHES 1,2 ET 3
!
    nedep = 0
    do j = 1, 3
        delta(j) = 1.d+10
        nefin = nedep+nbnozo(j)
        do inno = nedep+1, nefin
            noeu2 = tbnozo(inno)
            if (noeu2 .ne. noeu1) then
                dist = sqrt((coord(1)-xy(1, noeu2))**2+(coord(2)-xy(2, noeu2))**2)
                delta(j) = min(delta(j), dist)
            end if
        end do
        nedep = nefin
    end do
!
! 3 - CALCUL DE L ENERGIE POUR DIFFERENTS RAYONS RAYZ
!     ENER=SOMME(ENERGI(EF)*SPRIM(EF)/AIRE(EF))/AIRTOT
!     LA SOMME S EFFECTUE SUR TOUS LES EFS DES COUCHES 1,2 ET 3
!     AIRTOT EST L AIRE DU CERCLE DE CENTRE NOEU1 ET DE RAYON RAYZ
!     AIRE(EF) EST LA SURFACE DE L EF CONSIDERE
!     SPRIM(EF) EST LA SURFACE DE L EF CONSIDERE INCLUE DANS LE CERCLE
!     CAS 1 : EF INCLU DANS LE CERCLE => SPRIM(EF)=AIRE(EF)
!     CAS 2 : EF EXCLU DU CERCLE      => SPRIM(EF)=0
!     CAS 3 : EF INCLUE EN PARTIE     => SPRIM(EF) A CALCULER
!
    do iint = 1, nbint
!
        rayz(iint) = delta(1)+(delta(3)-delta(1))*(iint-1)/(nbint-1)
        ener(iint) = 0.d+0
        airtot = 0.d+0
!
! 3.1 - BOUCLE SUR LES EFS
!
        do inel = 1, nbelt
!
            nuef = tbelzo(inel)
            nsomm = icnc(1, nuef)
!
! 3.1.1 - NOMBRE DE NOEUD NINTER DE L EF INEL INCLU DANS LE CERCLE
!
            nint = 0
            ipoi1 = 0
            do inno = 1, nsomm
                coor(1, inno) = xy(1, icnc(inno+2, nuef))
                coor(2, inno) = xy(2, icnc(inno+2, nuef))
                ray = sqrt((coord(1)-coor(1, inno))**2+(coord(2)-coor(2, inno))**2)
                if (ray .le. rayz(iint)*epsir) then
                    nint = nint+1
                    if (ipoi1 .eq. 0) then
                        ipoi1 = inno
                    else
                        ipoi4 = inno
                    end if
                else
                    ipoi2 = inno
                end if
            end do
!
! 3.1.2 - AIRE DE L INTERSECTION SPRIM SELON LES CAS
!         SI NINT=NSOMM SPRIM = AIRE(EF)
!         SI NINT=0     SPRIM = 0
!         SINON ON CORRIGE
!
            sprim = 0.d+0
            if (nint .eq. nsomm) then
                sprim = aire(nuef)
!
! SI 2 NOEUDS APPARTIENNENT AU CERCLE SPRIM=AIRE - TRIANGLE COUPE
! RQ :DINTER CALCULE LES COORDONNEES DES NOEUDS COUPANT LE CERCLE
!
            else if (nint .eq. (nsomm-1)) then
                ip1 = ipoi2+1
                if (ip1 .gt. nsomm) ip1 = ip1-nsomm
                call dinter(coord, rayz(iint), coor(1, ipoi2), coor(1, ip1), coorin(1, 1))
                ip1 = ipoi2-1
                if (ip1 .le. 0) ip1 = ip1+nsomm
                call dinter(coord, rayz(iint), coor(1, ipoi2), coor(1, ip1), coorin(1, 2))
!
                call dcspri(coor(1, ipoi2), coorin, sprim)
                sprim = aire(nuef)-sprim
!
! SI 1 NOEUD APPARTIENT AU CERCLE SPRIM=TRIANGLE COUPE
!
            else if (nint .eq. 1) then
                ip1 = ipoi1+1
                if (ip1 .gt. 3) ip1 = ip1-3
                call dinter(coord, rayz(iint), coor(1, ipoi1), coor(1, ip1), coorin(1, 1))
                ip1 = ipoi1+2
                if (ip1 .gt. 3) ip1 = ip1-3
                call dinter(coord, rayz(iint), coor(1, ipoi1), coor(1, ip1), coorin(1, 2))
!
                call dcspri(coor(1, ipoi1), coorin, sprim)
!
! CAS PARTICULIER DES QUADRILATERES
!
            else if (nsomm .eq. 4 .and. nint .eq. 2) then
                if (ipoi1 .eq. 1 .and. ipoi4 .eq. 4) then
                    call dinter(coord, rayz(iint), coor(1, 1), coor(1, 2), coorin(1, 1))
                    call dinter(coord, rayz(iint), coor(1, 4), coor(1, 3), coorin(1, 2))
                    call dcqpri(coor(1, 1), coor(1, 4), coorin, sprim)
                else if (ipoi1 .eq. 3 .and. ipoi4 .eq. 4) then
                    call dinter(coord, rayz(iint), coor(1, 4), coor(1, 1), coorin(1, 1))
                    call dinter(coord, rayz(iint), coor(1, 3), coor(1, 2), coorin(1, 2))
                    call dcqpri(coor(1, 4), coor(1, 3), coorin, sprim)
                else if (ipoi1 .eq. 2 .and. ipoi4 .eq. 3) then
                    call dinter(coord, rayz(iint), coor(1, 3), coor(1, 4), coorin(1, 1))
                    call dinter(coord, rayz(iint), coor(1, 2), coor(1, 1), coorin(1, 2))
                    call dcqpri(coor(1, 3), coor(1, 2), coorin, sprim)
                else if (ipoi1 .eq. 1 .and. ipoi4 .eq. 2) then
                    call dinter(coord, rayz(iint), coor(1, 2), coor(1, 3), coorin(1, 1))
                    call dinter(coord, rayz(iint), coor(1, 1), coor(1, 4), coorin(1, 2))
                    call dcqpri(coor(1, 2), coor(1, 1), coorin, sprim)
                else
!             IPOI1 et IPOI4 non traité
                    ASSERT(.false.)
                end if
            end if
!
            ASSERT(sprim .ge. 0.d0)
!
! 3.1.3 - CALCUL DE L ENERGIE
!
            ener(iint) = ener(iint)+energi(nuef)*sprim/aire(nuef)
            airtot = airtot+sprim
!
! 2 - FIN DE LA BOUCLE SUR TOUS LES EF
!
        end do
!
        ASSERT(airtot .gt. 0.d0 .and. ener(iint) .gt. 0.d0)
        ener(iint) = ener(iint)/airtot
!
! LISSAGE DE LA COURBE ENERGI=F(RAYON) POUR IDENTIFIER PE
!
        if (iint .ge. 2) then
            i = iint-1
            ener(iint) = min(ener(iint), ener(i))
        end if
!
! 3 - FIN DE LA BOUCLE SUR LE CALCUL DE L ENERGIE
!
    end do
!
! 4 - CALCUL DU DEGRE DE LA SINGULARITE
!     PAR LA METHODE DES MOINDRES CARRES
!
    call dcalph(rayz, ener, nbint, pe)
!
end subroutine
