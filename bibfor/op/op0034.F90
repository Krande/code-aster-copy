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
subroutine op0034()
!
    implicit none
!
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/charth.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
! COMMAND:  AFFE_CHAR_THER_*
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
    character(len=3) :: exi_sech, exi_no_sech
    character(len=8) :: model
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Which command ?
    call getres(load, k16dummy, command)
    call getvid(' ', 'MODELE', scal=model)
    if (command .eq. 'AFFE_CHAR_THER') then
        valeType = 'REEL'
        call dismoi('EXI_SECH', model, 'MODELE', repk=exi_sech)
        if (exi_sech .eq. 'OUI') call utmess('F', 'CHARGES2_13')
    elseif (command .eq. 'AFFE_CHAR_SECH') then
        valeType = 'REEL'
        call dismoi('EXI_NON_SECH', model, 'MODELE', repk=exi_no_sech)
        if (exi_no_sech .eq. 'OUI') call utmess('F', 'CHARGES2_14')
    else if (command .eq. 'AFFE_CHAR_THER_F') then
        valeType = 'FONC'
        call dismoi('EXI_SECH', model, 'MODELE', repk=exi_sech)
        if (exi_sech .eq. 'OUI') call utmess('F', 'CHARGES2_13')
    else if (command .eq. 'AFFE_CHAR_SECH_F') then
        valeType = 'FONC'
        call dismoi('EXI_NON_SECH', model, 'MODELE', repk=exi_no_sech)
        if (exi_no_sech .eq. 'OUI') call utmess('F', 'CHARGES2_14')
    else
        ASSERT(ASTER_FALSE)
    end if

! - Load type
    call wkvect(load//'.TYPE', 'G V K8', 1, vk8=loadType)
    if (valeType .eq. 'REEL') then
        loadType(1) = 'THER_RE'
    else if (valeType .eq. 'FONC') then
        loadType(1) = 'THER_FO'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Loads treatment
    call charth(load, valeType)
!
    call jedema()
end subroutine
