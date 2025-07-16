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
subroutine ntload_chck(listLoad)
!
    use listLoad_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=24), intent(in) :: listLoad
!
! --------------------------------------------------------------------------------------------------
!
! THER_LINEAIRE - Algorithm
!
! Check loads
!
! --------------------------------------------------------------------------------------------------
!
! In  listLoad         : name of datastructure for list of loads
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, iLoad, nbLoad
    aster_logical :: isnotallowed
    character(len=8) :: loadName
    character(len=24) :: loadField
    character(len=24) :: loadNameJv
    character(len=24), pointer :: listLoadName(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    isnotallowed = .false.

! - Datastructure access
    call getNbLoadsFromList(listLoad, nbLoad)
    loadNameJv = listLoad(1:19)//'.LCHA'
    call jeveuo(loadNameJv, "E", vk24=listLoadName)

! - Seek for special loads
    call jeexin(loadNameJv, iret)
    if (iret .ne. 0) then
        if (nbLoad .ne. 0) then
            do iLoad = 1, nbLoad
                loadName = listLoadName(iLoad) (1:8)
                loadField = loadName(1:8)//'.CHTH'//'.CONVE'
                call jeexin(loadField(1:19)//'.VALE', iret)
                if (iret .ne. 0) then
                    isnotallowed = .true.
                end if
            end do
        end if
    end if
!
! - Fatal error
!
    if (isnotallowed) then
        call utmess('F', 'THERNONLINE4_1')
    end if
!
end subroutine
