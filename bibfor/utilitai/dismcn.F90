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

subroutine dismcn(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(CHAM_NO)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
!
#include "asterfort/dismgd.h"
#include "asterfort/dismeq.h"
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
    character(len=8) :: nogd
    integer(kind=8) :: iarefe, i8
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
    else if (questl .eq. 'NOM_MAILLA') then
        call jeveuo(nomob//'.REFE', 'L', iarefe)
        call dismeq(questl, zk24(iarefe-1+2), repi, repk, ierd)
    else if (questl .eq. 'NB_DDLACT') then
        call jeveuo(nomob//'.REFE', 'L', iarefe)
        call dismeq(questl, zk24(iarefe-1+2), repi, repk, ierd)
    else if (questl .eq. 'TYPE_CHAMP') then
        repk = 'NOEU'
    else if (questl(1:7) .eq. 'NUM_GD ') then
        call jeveuo(nomob//'.REFE', 'L', iarefe)
        call dismeq(questl, zk24(iarefe-1+2), repi, repk, ierd)
    else if (questl(1:7) .eq. 'NOM_GD ') then
        call jeveuo(nomob//'.REFE', 'L', iarefe)
        call dismeq(questl, zk24(iarefe-1+2), repi, repk, ierd)
    else if (questl .eq. 'TYPE_SUPERVIS' .or. questl .eq. 'TYPE_SCA') then
        call jeveuo(nomob//'.REFE', 'L', iarefe)
        call dismeq("NOM_GD", zk24(iarefe-1+2), i8, nogd, ierd)
        if (questl .eq. 'TYPE_SUPERVIS') then
            repk = 'CHAM_NO_'//nogd
        else
            call dismgd(questl, nogd, repi, repk, ierd)
        end if
    else if (questl .eq. 'NUME_EQUA') then
        call jeveuo(nomob//'.REFE', 'L', iarefe)
        repk = zk24(iarefe+1)
    else
        ierd = 1
    end if
!
    repkz = repk
    call jedema()
end subroutine
