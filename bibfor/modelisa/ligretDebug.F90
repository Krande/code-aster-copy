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
subroutine ligretDebug(ligretZ)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
!
    character(len=*), intent(in) :: ligretZ
!
! --------------------------------------------------------------------------------------------------
!
! Debug LIGRET
!
! In  ligret           : name of LIGRET
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: ligret
    character(len=8) :: docu, ligrelLgrf(3)
    character(len=16), pointer :: ligrelPhen(:) => null()
    character(len=16), pointer :: ligrelMode(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    ligret = ligretZ

    call dismoi('NOM_MAILLA', ligret, "LIGREL", repk=ligrelLgrf(1))
    call dismoi('PARTITION', ligret, "LIGREL", repk=ligrelLgrf(2))
    call dismoi('TYPE_LAGR', ligret, "LIGREL", repk=ligrelLgrf(3))

    call jelira(ligret//'.LGRF', "DOCU", cval=docu)
    call jeveuo(ligret//'.MODE', 'E', vk16=ligrelMode)
    call jeveuo(ligret//'.PHEN', 'E', vk16=ligrelPhen)

    WRITE (6, *) "================================================================================="
    WRITE (6, *) "Finite Element Descriptor for contact (LIGRET)"
    WRITE (6, *) "================================================================================="
    WRITE (6, *)
    WRITE (6, *) "Mesh : ", ligrelLgrf(1)
    WRITE (6, *) "Model : ", ligrelLgrf(2)
    WRITE (6, *) "Lagr12 : ", ligrelLgrf(3)
    WRITE (6, *) "SD Domjoints : ", ligrelLgrf(3)
    WRITE (6, *) "Phenomenon (docu) : ", docu
    WRITE (6, *) "Phenomenon : ", ligrelPhen(1)
    WRITE (6, *) "Modelisation : ", ligrelMode(1)

!
end subroutine
