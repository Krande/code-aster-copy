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
subroutine comp_meta_prnt(hasTemper, comporMetaInfo)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    aster_logical, intent(in) :: hasTemper
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
! In  hasTemper         : flag for tempering
! In  comporMetaInfo    : name of object for information about internal variables and comportement
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: mapZoneNume, mapNbZone, nbElemZone
    integer(kind=8) :: iVari, nbVari, ntVari, nbPhase
    character(len=16) :: metaLaw, metaType
    integer(kind=8), pointer :: comporInfoInfo(:) => null()
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    character(len=16), pointer :: comporInfoVari(:) => null()
    character(len=16), pointer :: comporInfoRela(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Access to informations
    call jeveuo(comporMetaInfo(1:19)//'.INFO', 'L', vi=comporInfoInfo)
    ntVari = comporInfoInfo(4)
    if (ntVari .eq. 0) then
        goto 99
    end if

    if (hasTemper) then
        call utmess('I', 'METALLURGY1_2')
    else
        call utmess('I', 'METALLURGY1_1')
    end if
    mapNbZone = comporInfoInfo(2)
    call jeveuo(comporMetaInfo(1:19)//'.RELA', 'L', vk16=comporInfoRela)
    call jeveuo(comporMetaInfo(1:19)//'.ZONE', 'L', vi=comporInfoZone)

    do mapZoneNume = 1, mapNbZone
        nbElemZone = comporInfoZone(mapZoneNume)
        if (nbElemZone .ne. 0) then
! --------- Acces to list of name of internal variables
            call jeveuo(jexnum(comporMetaInfo(1:19)//'.VARI', mapZoneNume), &
                        'L', vk16=comporInfoVari)
            call jelira(jexnum(comporMetaInfo(1:19)//'.VARI', mapZoneNume), 'LONMAX', nbVari)

! --------- Get names of relation
            metaType = comporInfoRela(3*(mapZoneNume-1)+1)
            metaLaw = comporInfoRela(3*(mapZoneNume-1)+2)
            read (comporInfoRela(3*(mapZoneNume-1)+3), '(I16)') nbPhase

! --------- Print name of internal variables
            call utmess('I', 'METALLURGY1_4', si=nbElemZone)
            call utmess('I', 'METALLURGY1_5', sk=metaType)
            call utmess('I', 'METALLURGY1_6', sk=metaLaw)
            call utmess('I', 'METALLURGY1_8', si=nbPhase)
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
