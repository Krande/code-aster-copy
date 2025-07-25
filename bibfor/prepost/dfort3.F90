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
subroutine dfort3(nsommx, icnc, noeu1, noeu2, tbelzo, &
                  nbelt, tbnozo, nbnoe, xy, volume, &
                  energi, pe)
!
!********************************************************************
!                  BUT DE CETTE ROUTINE :                           *
!      CALCULER LE DEGRE DE LA SINGULARITE PE SUR L ARETE           *
!                 D EXTREMITE NOEU1 ET NOEU2                        *
! 1) CALCUL DE L ENERGIE MOYENNE ENER SUR UN CYLINDRE D EXTREMITE   *
!    NOEU1 ET NOEU2 ET POUR DIFFERENTS RAYONS RAYZ                  *
! 2) COMPARAISON DE CETTE ENERGIE CALCULEE AVEC L ENERGIE THEORIQUE *
!    EN POINTE DE FISSURE ENER_THEO=K*(RAYZ**(2*(PE-1)))+C          *
!    ET CALCUL DU COEFFICENT PE PAR LA METHODE DES MOINDRES CARRES  *
!********************************************************************
!
! IN  NSOMMX               : NOMBRE DE SOMMETS MAX PAR EF
! IN  ICNC(NSOMMX+2,NELEM) : CONNECTIVITE EF => NOEUDS CONNECTES
!     1ERE VALEUR = NBRE DE NOEUDS SOMMETS CONNECTES A L EF N°X
!     2EME VALEUR = 1 SI EF UTILE 0 SINON
!     CONNECTIVITE  EF N°X=>N° DE NOEUDS SOMMETS CONNECTES A X
!     EN 2D EF UTILE = QUAD OU TRIA
!     EN 3D EF UTILE = TETRA OU HEXA
! IN  NOEU1 ET NOEU2       : NOEUDS DE L ARETE CONSIDEREE
! IN  TBELZO(NBELT)        : EFS DES COUCHES 1, 2 ET 3
! IN  NBELT                : NBR D EFS DES COUCHES 1,2,3
! IN  TBNOZO(NBNOE)        : NOEUDS DES COUCHES 1, 2 ET 3
! IN  NBNOE                : SOMME(NBNOZO)
! IN  XY(3,NNOEM)          : COORDONNEES DES NOEUDS
! IN  VOLUME(NELEM)        : VOLUME DE CHAQUE EF
! IN  ENERGI(NELEM)        : ENERGIE SUR CHAQUE EF
! OUT PE                   : DEGRE DE LA SINGULARITE
!
    implicit none
!
! DECLARATION GLOBALE
!
#include "asterfort/assert.h"
#include "asterfort/dcalph.h"
#include "asterfort/dinttc.h"
#include "asterfort/drao12.h"
#include "asterfort/dvolu1.h"
#include "asterfort/dvolu2.h"
#include "asterfort/dvolu3.h"
#include "asterfort/dvolu4.h"
    integer(kind=8) :: nsommx, icnc(nsommx+2, *), noeu1, noeu2
    integer(kind=8) :: nbelt, nbnoe
    integer(kind=8) :: tbelzo(nbelt), tbnozo(nbnoe)
    real(kind=8) :: xy(3, *), volume(*), energi(*), pe
!
! DECLARATION LOCALE
!
    integer(kind=8) :: i, inno, nuef, iint, inel, noeud
    integer(kind=8) :: nbint
    parameter(nbint=6)
    integer(kind=8) :: nint, nhop, npir, next, nbi, norm(2, 4)
    real(kind=8) :: coord1(3), coord2(3), coorn(3, 12)
    real(kind=8) :: xao1, yao1, zao1, xo1o2, yo1o2, zo1o2
    real(kind=8) :: do1o2
    real(kind=8) :: rayon, dist, ray(2), r, rayz(nbint)
    real(kind=8) :: voli, voltot
    real(kind=8) :: ener(nbint)
!
! 1 - COORDONNEES DE NOEU1 ET NOEU2
!
    do i = 1, 3
        coord1(i) = xy(i, noeu1)
        coord2(i) = xy(i, noeu2)
    end do
    xo1o2 = coord2(1)-coord1(1)
    yo1o2 = coord2(2)-coord1(2)
    zo1o2 = coord2(3)-coord1(3)
    do1o2 = sqrt(xo1o2**2+yo1o2**2+zo1o2**2)
    ASSERT(do1o2 .ge. 1.d-20)
!
! 2 - CALCUL DU RAYON MAXI
!
    rayon = -1.d+10
!
    do inno = 1, nbnoe
!
        noeud = tbnozo(inno)
        xao1 = xy(1, noeud)-coord1(1)
        yao1 = xy(2, noeud)-coord1(2)
        zao1 = xy(3, noeud)-coord1(3)
!
        dist = sqrt( &
               (yao1*zo1o2-zao1*yo1o2)**2+(zao1*xo1o2-xao1*zo1o2)**2+(xao1*yo1o2-yao1&
                                                                     &*xo1o2)**2 &
               )/do1o2
        rayon = max(rayon, dist)
!
    end do
!
! 3 - CALCUL DE L ENERGIE POUR DIFFERENTS RAYONS RAYZ
!     ENER=SOMME(ENERGI(EF)*VOLI/VOLUME(EF))/VOLTOT
!     LA SOMME S EFFECTUE SUR TOUS LES EFS DES COUCHES 1,2 ET 3
!     VOLTOT EST LE VOLUME DU CYLINDRE
!     VOLUME(EF) EST LE VOLUME DE L EF CONSIDERE
!     VOLI EST LE VOLUME DE L INTERSECTION ENTRE L EF CONSIDERE
!     ET LE CYLINDRE D EXTREMITE NOEU1 ET NOEU2 ET DE RAYON RAYZ
!
    do iint = 1, nbint
!
        r = (iint*rayon)/nbint
        rayz(iint) = r
        ener(iint) = 0.d+0
        voltot = 0.d+0
!
! 3.1 - BOUCLE SUR TOUS LES EFS
!
        do inel = 1, nbelt
!
            nuef = tbelzo(inel)
!
! 3.1.1 - CLASSEMENT DES NOEUDS
!         NINT NBR DE NOEUDS INTERNES AU CYLINDRE
!         NHOP NBR DE NOEUDS EXTERIEURS AUX DEUX PLANS
!         NPIR NBR DE NOEUDS INFEREURS AU RAYON
!         NEXT NBR DE NOEUDS EXTERNES AU CYLINDRE
!         NBI  NBR DE NOEUDS D INTERSECTION AU CYLINDRE
!
            nint = 0
            nhop = 0
            npir = 0
            next = 0
            do inno = 1, icnc(1, nuef)
                noeud = icnc(inno+2, nuef)
                call drao12(coord1, coord2, xo1o2, yo1o2, zo1o2, &
                            do1o2, xy(1, noeud), ray)
                if (ray(1) .le. rayz(iint)) then
                    if (ray(2) .eq. 0.0d0) then
                        norm(1, inno) = 1
                        norm(2, inno) = 0
                        nint = nint+1
                        npir = npir+1
                    else
                        norm(1, inno) = -1
                        if (ray(2) .eq. 1.0d0) norm(2, inno) = 1
                        if (ray(2) .eq. -1.0d0) norm(2, inno) = -1
                        next = next+10
                        nhop = nhop+1
                        npir = npir+1
                    end if
                else
                    if (ray(2) .ne. 0.0d0) then
                        norm(1, inno) = -1
                        if (ray(2) .eq. 1.0d0) norm(2, inno) = 1
                        if (ray(2) .eq. -1.0d0) norm(2, inno) = -1
                        next = next+1
                        nhop = nhop+1
                    else
                        norm(1, inno) = -1
                        norm(2, inno) = 0
                        next = next+1
                    end if
                end if
            end do
!
! 3.1.2 - VOLUME DE L INTERSECTION VOLI SELON LES CAS
!         SI NINT = 4, VOLI=VOLUME DE L EF
!         SI NEXT = 4, VOLI=0
!         SINON ON CORRIGE
!
            if (nint .eq. 4) then
                voli = volume(nuef)
            else if (next .ne. 4) then
                voli = 0.0d0
                do inno = 1, 4
                    coorn(1, inno) = xy(1, icnc(inno+2, nuef))
                    coorn(2, inno) = xy(2, icnc(inno+2, nuef))
                    coorn(3, inno) = xy(3, icnc(inno+2, nuef))
                end do
!
! CALCUL DES INTERSECTIONS DES ARETES DU TETRA AVEC LE CYLINDRE
!
                call dinttc(coord1, coord2, xo1o2, yo1o2, zo1o2, &
                            do1o2, r, norm, nint, nhop, &
                            npir, coorn, nbi)
!
                if (nint .eq. 1 .and. nbi .eq. 3) then
                    i = -1
                    voli = dvolu1(i, coorn, norm, volume)
                else if (nint .eq. 3 .and. nbi .eq. 3) then
                    voli = dvolu1(nuef, coorn, norm, volume)
                else if (nint .eq. 2 .and. nbi .eq. 2) then
                    voli = dvolu4(coorn, norm, coord1)
                else if (nint .eq. 2 .and. nbi .eq. 3) then
                    voli = dvolu3(coorn, norm, coord1)
                else if (nint .eq. 2 .and. nbi .eq. 4) then
                    voli = dvolu2(coorn, norm)
                else if (nbi .eq. 0) then
                    voli = 0.0d0
                else
!             problème intersection cylindre tétraèdre
                    ASSERT(.false.)
                end if
            end if
!
            ASSERT(voli .ge. 0.d0)
!
! 3.1.3 - CALCUL DE L ENERGIE
!
            ener(iint) = ener(iint)+energi(nuef)*voli/volume(nuef)
            voltot = voltot+voli
!
! 2 - FIN DE LA BOUCLE SUR TOUS LES EF
!
        end do
!
        ASSERT(voltot .gt. 0.d0 .and. ener(iint) .gt. 0.d0)
        ener(iint) = ener(iint)/voltot
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
!
! 4 - CALCUL DU DEGRE DE LA SINGULARITE
!     PAR LA METHODE DES MOINDRES CARRES
!
    call dcalph(rayz, ener, nbint, pe)
!
end subroutine
