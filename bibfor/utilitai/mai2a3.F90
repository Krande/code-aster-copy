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

subroutine mai2a3(mailla)
    implicit none

!METTRE UN MAILLAGE 2D EN 3D SI BESOIN APRES UN CHANGEMENT
!

#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"

    character(len=*) :: mailla
! IN  : MAILLA  : NOM DE LA SD MAILLAGE
    character(len=16) :: repk
    integer(kind=8), pointer :: dime(:) => null()

    call jemarq()

    call jeveuo(mailla//'.DIME', 'E', vi=dime)
    call dismoi('Z_ZERO', mailla, 'MAILLAGE', repk=repk)
    if (dime(6) == 2 .and. repk == 'NON') then
        dime(6) = 3
    end if

    call jedema()
end subroutine
