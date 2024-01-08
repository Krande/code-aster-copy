! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine comp_meta_prnt(comporMetaInfo)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(in) :: comporMetaInfo
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Print informations about internal variables
!
! --------------------------------------------------------------------------------------------------
!
! In  comporMetaInfo    : name of object for information about internal variables and comportement
!
! --------------------------------------------------------------------------------------------------
!
    integer :: iVari, mapZoneNume
    integer :: nbVari, mapNbZone, nb_elem_zone, nt_vari
    character(len=16) :: metaLaw, metaType
    integer, pointer :: comporInfoInfo(:) => null()
    integer, pointer :: comporInfoZone(:) => null()
    character(len=16), pointer :: comporInfoVari(:) => null()
    character(len=16), pointer :: comporInfoRela(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to informations
    call jeveuo(comporMetaInfo(1:19)//'.INFO', 'L', vi=comporInfoInfo)
    nt_vari = comporInfoInfo(4)
    if (nt_vari .eq. 0) then
        goto 99
    end if
    call utmess('I', 'METALLURGY1_1')
    mapNbZone = comporInfoInfo(2)
    call jeveuo(comporMetaInfo(1:19)//'.RELA', 'L', vk16=comporInfoRela)
    call jeveuo(comporMetaInfo(1:19)//'.ZONE', 'L', vi=comporInfoZone)

    do mapZoneNume = 1, mapNbZone
        nb_elem_zone = comporInfoZone(mapZoneNume)
        if (nb_elem_zone .ne. 0) then
! --------- Acces to list of name of internal variables
            call jeveuo(jexnum(comporMetaInfo(1:19)//'.VARI', mapZoneNume), &
                        'L', vk16=comporInfoVari)
            call jelira(jexnum(comporMetaInfo(1:19)//'.VARI', mapZoneNume), 'LONMAX', nbVari)

! --------- Get names of relation
            metaType = comporInfoRela(2*(mapZoneNume-1)+1)
            metaLaw = comporInfoRela(2*(mapZoneNume-1)+2)

! --------- Print name of internal variables
            call utmess('I', 'METALLURGY1_4', si=nb_elem_zone)
            call utmess('I', 'METALLURGY1_5', sk=metaType)
            call utmess('I', 'METALLURGY1_6', sk=metaLaw)
            call utmess('I', 'METALLURGY1_9', si=nbVari)
            do iVari = 1, nbVari
                call utmess('I', 'COMPOR4_20', sk=comporInfoVari(iVari), si=iVari)
            end do
        end if
    end do
!
99  continue
!
    call jedema()
!
end subroutine
