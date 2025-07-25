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
subroutine dismqu(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(DEGRE ELEMENTS 3D)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi
    character(len=32) :: repk
    character(len=19) :: nomob
    character(len=*) :: repkz, nomobz
! ----------------------------------------------------------------------
!    IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE LIGREL
!    OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=16) :: nomte
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iagrel, igr, iret, ite, n1, nbgr, nbnon
    integer(kind=8) :: nboui
!-----------------------------------------------------------------------
    call jemarq()
    nomob = nomobz
    repk = ' '
    repi = 0
    ierd = 0
!
    if (questi .eq. 'ELEM_VOLU_QUAD') then
        nboui = 0
        nbnon = 0
        call jeexin(nomob//'.LIEL', iret)
        ierd = 1
        if (iret .gt. 0) then
            call jelira(nomob//'.LIEL', 'NUTIOC', nbgr)
            do igr = 1, nbgr
                call jeveuo(jexnum(nomob//'.LIEL', igr), 'L', iagrel)
                call jelira(jexnum(nomob//'.LIEL', igr), 'LONMAX', n1)
                ite = zi(iagrel-1+n1)
                call jenuno(jexnum('&CATA.TE.NOMTE', ite), nomte)
                if (nomte .eq. 'MECA_HEXA20' .or. nomte .eq. 'MECA_HEXA27' .or. nomte .eq. &
                    'MECA_PENTA15' .or. nomte .eq. 'MECA_TETRA10' .or. nomte .eq. &
                    'MECA_PYRAM13' .or. nomte .eq. 'MECA_HEXS20' .or. nomte .eq. &
                    'MECA_PENTA18') then
                    repk = 'OUI'
                    nboui = nboui+1
                else if (nomte .eq. 'MECA_HEXA8' .or.&
 &          nomte .eq. 'MECA_PENTA6' .or. nomte .eq. 'MECA_TETRA4' .or.&
 &          nomte .eq. 'MECA_PYRAM5') then
                    repk = 'NON'
                    nbnon = nbnon+1
                end if
                ierd = 0
            end do
        end if
        if (nboui .ne. 0 .and. nbnon .ne. 0) repk = 'MEL'
    else
        ierd = 1
    end if
!
!
    repkz = repk
    call jedema()
end subroutine
