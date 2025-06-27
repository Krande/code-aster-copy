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
subroutine op0007()
!
    implicit none
!
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/charme.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! COMMAND:  AFFE_CHAR_MECA_*
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4) :: valeType
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
    if (command .eq. 'AFFE_CHAR_MECA') then
        valeType = 'REEL'
    else if (command .eq. 'AFFE_CHAR_MECA_F') then
        valeType = 'FONC'
    else if (command .eq. 'AFFE_CHAR_MECA_C') then
        valeType = 'COMP'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Load type
    call wkvect(load//'.TYPE', 'G V K8', 1, vk8=loadType)
    if (valeType .eq. 'REEL') then
        loadType(1) = 'MECA_RE'
    else if (valeType .eq. 'FONC') then
        loadType(1) = 'MECA_FO'
    else if (valeType .eq. 'COMP') then
        loadType(1) = 'MECA_RI'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Loads treatment
    call charme(load, valeType)
!
    call jedema()
end subroutine
