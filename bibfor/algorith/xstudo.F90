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

subroutine xstudo(ndime, ninter, npts, nptm, ainter, &
                  nbpi, ip1, ip2, pm1a, pm1b, pm2)
    implicit none
!
#    include "jeveux.h"
#    include "asterfort/assert.h"
#    include "asterfort/xxmmvd.h"
    integer(kind=8) :: ndime, ninter, npts, nptm
    real(kind=8) :: ainter(*)
    integer(kind=8) :: nbpi, ip1(4), ip2(4), pm1a(4), pm1b(4), pm2(4)
!    BUT: CONSTRUIRE LES FACES DE DECOUPE DU TETRA DE DECOUPE DE REFERENCE
!         CALCULE LES POINTS D INTERSECTION ET POINTS MILIEUX PRESENTS PAR FACE
!
!     ENTREE
!       NINTER  : NOMBRE DE POINTS D'INTERSECTION
!       NPTS    : NOMBRE DE POINTS D'INTERSECTION CONFONDU AVEC UN NOEUD SOMMET
!       NPTM    : NOMBRE DE POINTS D'INTERSECTION CONFONDU AVEC UN NOEUD MILIEU
!       FPM1    : NUMRO DU POINT D'INTERSECTION CORRESPONDANT AU NOEUD MILIEU A NE PAS CONSIDERER
!                 DANS LE CAS NINTER=4, NPTS=2, NPTM=1
!     SORTIE
!       NBPI    : NOMBRE DE FACE INTERSECTE
!       IP1     : PREMIER IP DANS DANS UN FACE
!       IP2     : DEUXIEME IP DANS UN FACE
!       PM1A    : PREMIER PM1 DANS UN FACE
!       PM1B    : DEUXIEME PM1 DANS UN FACE
!       PM2     : PM2 DANS UN FACE
!     ----------------------------------------------------------------
!
!
! --------------------------------------------------------------------
!
!   ATTENTION:: (IP1,PM1A) (IP2,PM1B) NE DOIVENT PAS ETRE INTERVERTIS
!               SI LE CALCUL DE XMILFA CALCULE LE MILIEU D UNE ARETE DE DECOUPE
!               AU LIEU D APPELER XNEWTO POUR RETROUVER LE CENTRE DU QUAD8
!               L ORDRE PROPOSE ICI DEVIENT CORRELE AU SOUS DECOUPAGE XPENTE DANS XDECQV
!
    integer(kind=8) :: zxain

    zxain = xxmmvd('ZXAIN')

    if (ninter .eq. 1) then
        nbpi = 0
    else if (ninter .eq. 2 .and. ndime .eq. 2) then
        if (npts .eq. 0) then
            nbpi = 1
!       MILIEU DE I1-I2
            ip1(1) = 1
            ip2(1) = 2
            pm1a(1) = 1
            pm1b(1) = 3
            pm2(1) = 5
        else if (npts .eq. 1) then
            nbpi = 1
!       MILIEU DE I1-I2
            ip1(1) = 1
            ip2(1) = 2
        end if
    else if (ninter .eq. 3 .and. ndime .eq. 2) then
        nbpi = 1
!       MILIEU DE I1-I2
        ip1(1) = 2
        ip2(1) = 3
        pm1a(1) = 1
        pm1b(1) = 3
        pm2(1) = 5
    else if (ninter .eq. 3 .and. ndime .eq. 3) then
        if (npts .eq. 0) then
            nbpi = 3
!       MILIEU DE I1-I2
            ip1(1) = 1
            ip2(1) = 2
            pm1a(1) = 1
            pm1b(1) = 3
            pm2(1) = 7
!       MILIEU DE I2-I3
            ip1(2) = 2
            ip2(2) = 3
            pm1a(2) = 3
            pm1b(2) = 5
            pm2(2) = 8
!       MILIEU DE I3-I1
            ip1(3) = 1
            ip2(3) = 3
            pm1a(3) = 1
            pm1b(3) = 5
            pm2(3) = 9
        else if (npts .eq. 1) then
            nbpi = 3
!       MILIEU DE I1-I2
            ip1(1) = 1
            ip2(1) = 2
            pm1a(1) = 0
            pm1b(1) = 0
            pm2(1) = 5
!       MILIEU DE I2-I3
            ip1(2) = 2
            ip2(2) = 3
            pm1a(2) = 1
            pm1b(2) = 3
            pm2(2) = 6
!       MILIEU DE I3-I1
            ip1(3) = 3
            ip2(3) = 1
            pm1a(3) = 0
            pm1b(3) = 0
            pm2(3) = 7
        end if
    else if (ninter .eq. 4 .and. ndime .eq. 3) then
!       CAS 1 : DECOUPAGE DU TETRA EN DEUX PENTAEDRES
        if (npts .eq. 0) then
            nbpi = 4
            !       MILIEU DE I1-I2
            ip1(1) = 2
            ip2(1) = 1
            pm1a(1) = 3
            pm1b(1) = 1
            pm2(1) = 9
            !       MILIEU DE I2-I4
            ip1(2) = 4
            ip2(2) = 2
            pm1a(2) = 8
            pm1b(2) = 4
            pm2(2) = 10
            !       MILIEU DE I4-I3
            ip1(3) = 4
            ip2(3) = 3
            pm1a(3) = 7
            pm1b(3) = 5
            pm2(3) = 11
            !       MILIEU DE I3-I1
            ip1(4) = 3
            ip2(4) = 1
            pm1a(4) = 6
            pm1b(4) = 2
            pm2(4) = 12
!       CAS 2 : DECOUPAGE DU TETRA EN DEUX PENTAEDRES
        else if (npts .eq. 2 .and. nptm .ge. 1) then
!    LE POINT D INTERSECTION NON CONFONDU AVEC UN NOEUD SOMMET EST EN 3EME POSITION
!    VOIRE XDECQU
            ASSERT(nint(ainter(zxain*(3-1)+1)) .ne. 0)
            nbpi = 2
            !       MILIEU DE I1-P
            ip1(1) = 1
            ip2(1) = 3
            !       MILIEU DE I2-P
            ip1(2) = 2
            ip2(2) = 3
        else if (npts .eq. 1) then
            ASSERT(nptm .ge. 1)
            nbpi = 3
!         MILIEU DE I1-I2
            ip1(1) = 2
            ip2(1) = 3
            pm2(1) = 7
!         MILIEU DE I2-I3
            ip1(2) = 3
            ip2(2) = 4
            pm2(2) = 8
!         MILIEU DE I3-I1
            ip1(3) = 4
            ip2(3) = 1
        else
            ASSERT(.false.)
        end if
    else
        ASSERT(.false.)
!
    end if
end subroutine
