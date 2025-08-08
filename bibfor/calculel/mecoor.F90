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

subroutine mecoor(ligrel, chgeom)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    character(len=19) :: ligrel
    character(len=19) :: chgeom
!
! ----------------------------------------------------------------------
!
! RETOURNE LE CHAMP DES COORDONNEES DU MAILLAGE
!
! ----------------------------------------------------------------------
!
!
! IN  NOMO   : NOM DU LIGREL
! OUT CHGEOM : CHAMP DE GEOMETRIE
!
!
!
!
    character(len=8) :: ma
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=ma)
    chgeom = ma//'.COORDO'
!
    call jedema()
!
end subroutine
