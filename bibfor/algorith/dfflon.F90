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

subroutine dfflon(geom, nonoff, nomnoe, inoff, nbnoff, &
                  typfon, d)
!
    implicit none
!
#include "jeveux.h"
!
#include "asterfort/dis2no.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/char8_to_int.h"
    real(kind=8) :: geom(*), d
    integer(kind=8) :: inoff, nbnoff
    character(len=8) :: nonoff(*), typfon
    character(len=24) :: nomnoe
!
!
! FONCTION REALISEE (OPERATEURS DEFI_FOND_FISS) :
!
!      RETOURNE UNE ESTIMATION DE LA LONGUEUR DES SEGMENTS DU FOND DE
!      FISSURE EN UN NOEUD DU FOND
!      (UNIQUEMENT EN 3D)
!
! IN
!   GEOM    : VALEURS DES COORDONNEES DU MAILLAGE
!   NONOFO  : NOMS DES NOEUDS DU FOND
!   NOMNOE  : OBJET '.NOMNOE' DU MAILLAGE
!   INOFF   : INCIDE LOCAL DU NOEUD DE LA BASE DEMANDE
!   NBNOFF  : NOMBRE DE NOEUDS DU FOND DE FISSURE
!   TYPFON  : TYPE DE FOND (OUVERT/FERME)
!
! OUT
!   D   : LONGUEUR CARACTERISTIQUES DES SEGMENTS DU FOND AUTOUT DU NOEUD
!
!-----------------------------------------------------------------------
!
!
    real(kind=8) :: dij, dih
    integer(kind=8) :: nunoi, nunoj, nunoh
!
!-----------------------------------------------------------------------
!
    call jemarq()
!
!     REMARQUE GENERALE : ON POURRAIT TROUVER UN ALGO PLUS ELEGANT,
!                         MAIS ON PERDRAIT EN LISIBILITE
!
!     NUNOI : NUMERO DU NOEUD COURANT
!     NUNOJ : NUMERO DU NOEUD APRES
!     NUNOH : NUMERO DU NOEUD AVANT
!
!     CAS PARTICULIER DU PREMIER NOEUD DU FOND DE FISSURE
    if (inoff .eq. 1) then
!
        nunoi = char8_to_int(nonoff(inoff))
        nunoj = char8_to_int(nonoff(inoff+1))
!
        dij = dis2no(geom, nunoi, nunoj)
        d = dij
!
        if (typfon .eq. 'FERME') then
!         ATTENTION, LE DERNIER NOEUD (EN POSITION NBNOFF) = LE 1ER
            nunoi = char8_to_int(nonoff(nbnoff))
            nunoh = char8_to_int(nonoff(nbnoff-1))
            dih = dis2no(geom, nunoi, nunoh)
            d = min(dij, dih)
        end if
!
!     CAS PARTICULIER DU DENIER NOEUD DU FOND DE FISSURE
    else if (inoff .eq. nbnoff) then
!
        nunoi = char8_to_int(nonoff(inoff))
        nunoh = char8_to_int(nonoff(inoff-1))
!
        dih = dis2no(geom, nunoi, nunoh)
        d = dih
!
        if (typfon .eq. 'FERME') then
!         ATTENTION, LE PREMIER NOEUD = LE DERNIER
            nunoi = char8_to_int(nonoff(1))
            nunoj = char8_to_int(nonoff(2))
            dij = dis2no(geom, nunoi, nunoj)
            d = min(dij, dih)
        end if
!
!     CAS GENERAL
    else
!
        nunoh = char8_to_int(nonoff(inoff-1))
        nunoi = char8_to_int(nonoff(inoff))
        nunoj = char8_to_int(nonoff(inoff+1))
!
        dih = dis2no(geom, nunoi, nunoh)
        dij = dis2no(geom, nunoi, nunoj)
        d = min(dij, dih)
!
    end if
!
    call jedema()
end subroutine
