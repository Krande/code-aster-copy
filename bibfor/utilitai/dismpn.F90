! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine dismpn(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(PROF_CHNO)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismlg.h"
#include "asterfort/dismcn.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/jexnom.h"
!
    integer :: repi, ierd
    character(len=*) :: questi
    character(len=*) :: nomobz, repkz
    character(len=32) :: repk
    character(len=19) :: nomob
! ----------------------------------------------------------------------
!     IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE PROF_CHNO
!     OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
!     LISTE DES QUESTIONS ADMISSIBLES:
!        'NB_DDLACT'
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=19) :: noligr
!
!
!
!-----------------------------------------------------------------------
    integer :: i, nbddlb, nbnos, nequ, nlili
!-----------------------------------------------------------------------
    call jemarq()
    nomob = nomobz
    repk = ' '
    repi = 0
    ierd = 0
!
    if (questi .eq. 'NB_DDLACT') then
!     --------------------------------
        call jelira(nomob//'.NUEQ', 'LONMAX', nequ)
        call jelira(nomob//'.LILI', 'NUTIOC', nlili)
        nbddlb = 0
        do i = 2, nlili
            call jenuno(jexnum(nomob//'.LILI', i), noligr)
            call dismlg('NB_NO_SUP', noligr, nbnos, repk, ierd)
            nbddlb = nbddlb+nbnos
        end do
        repi = nequ-3*(nbddlb/2)
!
!
    else if (questi .eq. 'NB_EQUA') then
!     --------------------------------
        call jelira(nomob//'.NUEQ', 'LONMAX', repi)
!
    else
        ierd = 1
    end if
!
    repkz = repk
    call jedema()
end subroutine
