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
subroutine load_unde_diri(listLoadZ)
!
    use listLoad_module
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/copich.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedup1.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: listLoadZ
!
! --------------------------------------------------------------------------------------------------
!
! List of loads - Utility
!
! Special copy for undead Dirichlet loads
!
! --------------------------------------------------------------------------------------------------
!
! In  listLoad          : name of datastructure for list of loads
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbLoad, loadNumeDiri, iLoad, iLoadDiri, i
    character(len=8) :: loadName, load_name_loca
    character(len=24) :: loadInfoJv, loadNameJv
    integer(kind=8), pointer :: listLoadInfo(:) => null()
    character(len=24), pointer :: listLoadName(:) => null()
    character(len=24), pointer :: noli(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Datastructure access
    loadInfoJv = listLoadZ(1:19)//'.INFC'
    loadNameJv = listLoadZ(1:19)//'.LCHA'
    call jeveuo(loadNameJv, 'E', vk24=listLoadName)
    call jeveuo(loadInfoJv, 'E', vi=listLoadInfo)
!
    call getNbLoadsFromList(listLoadZ, nbLoad)
    ASSERT(nbLoad .gt. 0)
    iLoadDiri = 0

! - Copy loads
    do iLoad = 1, nbLoad
        loadNumeDiri = listLoadInfo(iLoad+1)
        if (loadNumeDiri .eq. 4) then
            iLoadDiri = iLoadDiri+1
            load_name_loca = '&&DIRISU'
            if (iLoadDiri .gt. 99) then
                call utmess('F', 'CHARGES_29')
            end if
            call codent(iLoadDiri, 'D0', load_name_loca(7:8))
            loadName = listLoadName(iLoad) (1:8)

!           -- delete :
            call jedetr(load_name_loca//'.TYPE')
            call jedetr(load_name_loca//'.CHME.MODEL.NOMO')
            call detrsd('LIGREL', load_name_loca//'.CHME.LIGRE')
            call detrsd('CARTE', load_name_loca//'.CHME.CMULT')
            call detrsd('CARTE', load_name_loca//'.CHME.CIMPO')
            call jedetr(load_name_loca//'.DUAL.PRDK')
            call jedetr(load_name_loca//'.DUAL.PRDI')
            call jedetr(load_name_loca//'.DUAL.PRDSO')

!           -- copy :
            call jedup1(loadName//'.TYPE', 'V', load_name_loca//'.TYPE')
            call jedup1(loadName//'.CHME.MODEL.NOMO', 'V', load_name_loca//'.CHME.MODEL.NOMO')
            call copisd('LIGREL', 'V', loadName//'.CHME.LIGRE', load_name_loca//'.CHME.LIGRE')
            call copich('V', loadName//'.CHME.CMULT', load_name_loca//'.CHME.CMULT')
            call copich('V', loadName//'.CHME.CIMPO', load_name_loca//'.CHME.CIMPO')
            call jeveuo(load_name_loca//'.CHME.CIMPO.NOLI', 'E', vk24=noli)
            do i = 1, size(noli)
                if (noli(i) (1:8) .eq. loadName) then
                    noli(i) (1:8) = load_name_loca
                end if
            end do
            call jeveuo(load_name_loca//'.CHME.CMULT.NOLI', 'E', vk24=noli)
            do i = 1, size(noli)
                if (noli(i) (1:8) .eq. loadName) then
                    noli(i) (1:8) = load_name_loca
                end if
            end do
            call jedup1(loadName//'.DUAL.PRDK', 'V', load_name_loca//'.DUAL.PRDK')
            call jedup1(loadName//'.DUAL.PRDI', 'V', load_name_loca//'.DUAL.PRDI')
            call jedup1(loadName//'.DUAL.PRDSO', 'V', load_name_loca//'.DUAL.PRDSO')

            listLoadName(iLoad) (1:8) = load_name_loca
            listLoadName(iLoad) (9:16) = loadName
        end if
    end do
!
    call jedema()
!
end subroutine
