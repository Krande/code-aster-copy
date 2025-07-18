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

subroutine cfconn(defico, jdecno, ino, posno)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=24) :: defico
    integer(kind=8) :: jdecno, ino, posno
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES MAILLEES - UTILITAIRE)
!
! NOEUD ATTACHE A LA MAILLE (CONNECTIVITE DIRECTE)
!
! ----------------------------------------------------------------------
!
!
! IN  DEFICO : SD DE CONTACT (DEFINITION)
! IN  INO    : NUMERO ORDRE DU NOEUD DANS SD CONN.
! IN  JDECNO : DECALAGE POUR LECTURE DANS SD CONN.
! OUT POSNO  : POSITION DU NOEUD
!
!
!
!
    character(len=24) :: nomaco
    integer(kind=8) :: jnoma
    integer(kind=8) :: numno
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- ACCES SD DE CONTACT
!
    nomaco = defico(1:16)//'.NOMACO'
    call jeveuo(nomaco, 'L', jnoma)
!
! --- REPONSE
!
    if (jdecno .eq. -1) then
        ASSERT(.false.)
    end if
    numno = jdecno+ino
    posno = zi(jnoma+numno-1)
!
    call jedema()
!
end subroutine
