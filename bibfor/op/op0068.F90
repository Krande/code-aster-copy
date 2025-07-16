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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine op0068()
!
    implicit none
!
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/charac.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! COMMAND:  AFFE_CHAR_ACOU
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: load
    character(len=16) :: command, k16dummy
    character(len=8), pointer :: loadType(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Which command ?
    call getres(load, k16dummy, command)
    ASSERT(command .eq. 'AFFE_CHAR_ACOU')

! - Load type
    call wkvect(load//'.TYPE', 'G V K8', 1, vk8=loadType)
    loadType(1) = 'ACOU_RE'

! - Loads treatment
    call charac(load)
!
    call jedema()
end subroutine
