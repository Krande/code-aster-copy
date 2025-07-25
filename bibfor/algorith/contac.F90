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

subroutine contac(macor, nbcor, macoc, nbcoc, lface,&
                  lomodi, locorr, loreor, ma)
!
!  ROUTINE CONTAC
!    ROUTINE D ORIENTATION EN FONCTION DE KTYC
!  DECLARATIONS
!    KTYC   : NOM DU TYPE                   POUR UNE MAILLE FISSURE
!    KTYR   : NOM DU TYPE                   POUR UNE MAILLE REFERENCE
!    LOMODI : LOGICAL PRECISANT SI LA MAILLE EST UNE MAILLE MODIFIE
!    LOCORR : LOGICAL PRECISANT SI LA MAILLE EST BIEN ORIENTEE
!    LOREOR : LOGICAL PRECISANT SI LA MAILLE EST REORIENTEE
!    MA     : L OBJET DU MAILLAGE
!    MACOC  : TABLEAU DES NOMS DES NOEUDS   POUR UNE MAILLE FISSURE
!    MACOR  : TABLEAU DES NOMS DES NOEUDS   POUR UNE MAILLE REFERENCE
!    NBCOC  : NOMBRE DE CONNEX              POUR UNE MAILLE FISSURE
!    NBCOR  : NOMBRE DE CONNEX              POUR UNE MAILLE REFERENCE
!
!  MOT_CLEF : ORIE_FISSURE
!
!
    implicit none
!
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "asterfort/conhex.h"
#include "asterfort/conpen.h"
#include "asterfort/conqua.h"
    character(len=8) :: ktyc, ktyr
    integer(kind=8) :: nbcoc, nbcor
    character(len=8) :: macor(nbcor+2), macoc(nbcoc+2), ma
!
    aster_logical :: lface, lomodi, locorr, loreor
!
!-----------------------------------------------------------------------
!     FONCTIONS FORMULES PERMETTANT DE SAVOIR SI L'APPUI EST POSSIBLE
#define qua() (ktyc.eq.'QUAD4'.and.(ktyr.eq.'QUAD4'.or.ktyr.eq.'TRIA3')) \
    .or. (ktyc.eq.'QUAD8'.and.(ktyr.eq.'QUAD9'.or.ktyr.eq.'QUAD8' \
    .or.ktyr.eq.'TRIA6'))
#define pen() (ktyc.eq.'PENTA6 '.and. \
    (ktyr.eq.'PENTA6 '.or.ktyr.eq.'TETRA4')) \
    .or. (ktyc.eq.'PENTA15'.and. \
    (ktyr.eq.'PENTA15'.or.ktyr.eq.'TETRA10'))
#define hex() (ktyc.eq.'HEXA8 '.and. \
    (ktyr.eq.'HEXA8 '.or.ktyr.eq.'PENTA6 '.or.ktyr.eq.'PYRAM5 ')) \
    .or. (ktyc.eq.'HEXA20'.and. \
    (ktyr.eq.'HEXA20'.or.ktyr.eq.'PENTA15'.or.ktyr.eq.'PYRAM13'))
!     ------------------------------------------------------------------
!
    ktyc = macoc(2)
    ktyr = macor(2)
!
    if (qua()) then
!
        call conqua(macor, nbcor, macoc, nbcoc, lface,&
                    lomodi, locorr, loreor, ma)
!
    else if (pen()) then
!
        call conpen(macor, nbcor, macoc, nbcoc, lface,&
                    locorr, loreor, ma)
!
    else if (hex()) then
!
        call conhex(macor, nbcor, macoc, nbcoc, lface,&
                    lomodi, locorr, loreor, ma)
!
    else
!
        loreor=.false.
!
    endif
!
!     ------------------------------------------------------------------
end subroutine
