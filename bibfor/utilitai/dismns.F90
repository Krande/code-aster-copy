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

subroutine dismns(questi, nomobz, repi, repkz, ierd)
    implicit none
#include "jeveux.h"
!
#include "asterfort/dismgd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi, nomobz, repkz
! person_in_charge: jacques.pellet at edf.fr
!
!     --     DISMOI(CHAM_NO_S)
!
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE LIGREL
!
! OUT : REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, -1 --> CHAMP INEXISTANT)
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: iret, gd
    character(len=8) :: nogd
    character(len=19) :: nomob
    character(len=24) :: questl
    character(len=32) :: repk
    character(len=8), pointer :: cnsk(:) => null()
! DEB-------------------------------------------------------------------
!
    call jemarq()
    repk = ' '
    repi = 0
    ierd = 0
!
!
    nomob = nomobz
    questl = questi
!
    call jeexin(nomob//'.CNSK', iret)
    if (iret .eq. 0) then
        ierd = 1
        goto 9999
    end if
!
    call jeveuo(nomob//'.CNSK', 'L', vk8=cnsk)
    nogd = cnsk(2)
    call jenonu(jexnom('&CATA.GD.NOMGD', nogd), gd)
!
    if (questi .eq. 'TYPE_CHAMP') then
        repk = 'CNOS'
!
    else if (questi .eq. 'NOM_MAILLA') then
        repk = cnsk(1)
!
    else if (questl(1:6) .eq. 'NUM_GD') then
        repi = gd
!
    else if (questl(1:6) .eq. 'NOM_GD') then
        repk = nogd
!
    else if (questi .eq. 'TYPE_SCA') then
        call dismgd(questi, nogd, repi, repk, ierd)
!
    else
        ierd = 1
    end if
!
9999 continue
    repkz = repk
!
    call jedema()
end subroutine
