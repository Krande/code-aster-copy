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
subroutine xstami(noma, nmafon, nmaen1, nmaen2, nmaen3, &
                  jmafon, jmaen1, jmaen2, jmaen3)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nmafon, nmaen1, nmaen2, nmaen3
    integer(kind=8) :: jmafon, jmaen1, jmaen2, jmaen3
    character(len=8) :: noma
!
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM
!
! IMPRESSION DU STATUT DES MAILLES
!
! ----------------------------------------------------------------------
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NMAFON : NOMBRE DE MAILLES CONTENANT LE FOND DE FISSURE
! IN  NMAEN1 : NOMBRE DE MAILLES 'HEAVISIDE'
! IN  NMAEN2 : NOMBRE DE MAILLES 'CRACKTIP'
! IN  NMAEN3 : NOMBRE DE MAILLES 'HEAVISIDE-CRACKTIP'
! IN  JMAFON : POINTEUR SUR MAILLES 'CONTENANT LE FOND DE FISSURE
! IN  JMAEN1 : POINTEUR SUR MAILLES 'HEAVISIDE'
! IN  JMAEN2 : POINTEUR SUR MAILLES 'CRACKTIP'
! IN  JMAEN3 : POINTEUR SUR MAILLES 'HEAVISIDE-CRACKTIP'
!
!
!
!
    integer(kind=8) :: ifm, niv, ima
    character(len=8) :: nomail
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
!
    if (niv .ge. 2) then
        call utmess('I', 'XFEM_29', si=nmafon)
    end if
    if (niv .ge. 3) then
        do ima = 1, nmafon
            nomail = int_to_char8(zi(jmafon-1+ima))
            write (ifm, *) 'MAILLE ', nomail
        end do
    end if
!
    if (niv .ge. 2) then
        call utmess('I', 'XFEM_30', si=nmaen1)
    end if
    if (niv .ge. 3) then
        do ima = 1, nmaen1
            nomail = int_to_char8(zi(jmaen1-1+ima))
            write (ifm, *) 'MAILLE ', nomail
        end do
    end if
!
    if (niv .ge. 2) then
        call utmess('I', 'XFEM_31', si=nmaen2)
    end if
    if (niv .ge. 3) then
        do ima = 1, nmaen2
            nomail = int_to_char8(zi(jmaen2-1+ima))
            write (ifm, *) 'MAILLE ', nomail
        end do
    end if
!
    if (niv .ge. 2) then
        call utmess('I', 'XFEM_32', si=nmaen3)
    end if
    if (niv .ge. 3) then
        do ima = 1, nmaen3
            nomail = int_to_char8(zi(jmaen3-1+ima))
            write (ifm, *) 'MAILLE ', nomail
        end do
    end if
!
    call jedema()
end subroutine
