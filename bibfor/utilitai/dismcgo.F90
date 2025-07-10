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

subroutine dismcgo(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(CHAM_GEOM)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
!
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi
    character(len=*) :: nomobz, repkz
    character(len=24) :: questl
    character(len=32) :: repk
    character(len=19) :: nomob
! ----------------------------------------------------------------------
!     IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOB  : NOM D'UN OBJET DE TYPE NUM_DDL
!     OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPK   : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    integer(kind=8) :: iadesc
!-----------------------------------------------------------------------
    call jemarq()
    repk = ' '
    repi = 0
    ierd = 0
!
    nomob = nomobz
    questl = questi
!
    if (questl .eq. 'NB_EQUA') then
        call jelira(nomob//'.VALE', 'LONMAX', repi)
    else if (questl .eq. 'TYPE_CHAMP') then
        repk = 'GEOM'
    else if (questl(1:7) .eq. 'NUM_GD ') then
        call jeveuo(nomob//'.DESC', 'L', iadesc)
        repi = zi(iadesc)
    else if (questl .eq. "NB_NO_MAILLA") then
        call jelira(nomob//'.VALE', 'LONMAX', repi)
        repi = repi/3
    else if (questl(1:7) .eq. 'NOM_GD ') then
        repk = "GEOM_R"
    else
        ierd = 1
    end if
!
    repkz = repk
    call jedema()
end subroutine
