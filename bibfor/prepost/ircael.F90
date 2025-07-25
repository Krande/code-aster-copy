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

subroutine ircael(jcesdi, jcesli, jcesvi, jcesci, nummai, &
                  nbqcou, nbtcou, nbrsec, nbrfib, nbrgrf, nugrfi)
    implicit none
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
    integer(kind=8) :: nmaxfi
    parameter(nmaxfi=10)
    integer(kind=8) :: jcesdi, jcesli, nummai, nbqcou, nbtcou, jcesvi
    integer(kind=8) :: nbrsec, nbrfib, nbrgrf, nugrfi(nmaxfi), jcesci
! person_in_charge: nicolas.sellenet at edf.fr
! ----------------------------------------------------------------------
!  IMPR_RESU - CARACTERISTIQUE DE L'ELEMENT
!  -    -      --                   --
! ----------------------------------------------------------------------
!
!  CETTE ROUTINE EST UTILE POUR TROUVER LE NOMBRE DE COUCHES, ...
!    PRESENTES SUR UNE MAILLE
!
! IN  :
!   JCESDI  I    ADRESSE JEVEUX DU CESD
!   JCESLI  I    ADRESSE JEVEUX DU CESL
!   JCESVI  I    ADRESSE JEVEUX DU CESV
!   JCESCI  I    ADRESSE JEVEUX DU CESC
!   NUMMAI  I    NUMERO DE MAILLE
!
! OUT :
!   NBRCOU  I    NOMBRE DE COUCHES
!   NBRSEC  I    NOMBRE DE SECTEURS
!   NBRFIB  I    NOMBRE DE FIBRES
!   NBRGRF  I    NOMBRE DE GROUPES DE FIBRES
!   NUGRFI  I    NUMERO DU GROUPE DE FIBRES
!
#include "jeveux.h"
!
    integer(kind=8) :: nbrcmp, numgrf, icmp, iad
!
    nugrfi(1:nmaxfi) = 0
!
    nbqcou = 0
    nbtcou = 0
    nbrsec = 0
    nbrfib = 0
    nbrgrf = 0
    nugrfi = 0
!
    nbrcmp = zi(jcesdi-1+5+4*(nummai-1)+3)
    numgrf = 1
    do icmp = 1, nbrcmp
        call cesexi('C', jcesdi, jcesli, nummai, 1, 1, icmp, iad)
!
        if (iad .gt. 0) then
            if (zk8(jcesci+icmp-1) .eq. 'COQ_NCOU') then
                nbqcou = zi(jcesvi-1+iad)
            else if (zk8(jcesci+icmp-1) .eq. 'TUY_NCOU') then
                nbtcou = zi(jcesvi-1+iad)
            else if (zk8(jcesci+icmp-1) .eq. 'TUY_NSEC') then
                nbrsec = zi(jcesvi-1+iad)
            else if (zk8(jcesci+icmp-1) .eq. 'NBFIBR') then
                nbrfib = zi(jcesvi-1+iad)
            else if (zk8(jcesci+icmp-1) .eq. 'NBGRFI') then
                nbrgrf = zi(jcesvi-1+iad)
            else if (zk8(jcesci+icmp-1) (1:3) .eq. 'NUG') then
                nugrfi(numgrf) = zi(jcesvi-1+iad)
                numgrf = numgrf+1
            end if
        end if
    end do
    if (nbrfib .ne. 0) then
        ASSERT(nbrgrf .ne. 0)
    end if
!
end subroutine
