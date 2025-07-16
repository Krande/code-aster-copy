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
subroutine acdoch(listLoadPrep, listLoadZ, jvBase)
!
    use listLoad_module
    use listLoad_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "LoadTypes_type.h"
!
    type(ListLoad_Prep), intent(in) :: listLoadPrep
    character(len=*), intent(in) :: listLoadZ
    character(len=1), intent(in) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
! Acoustic - Read parameters
!
! Get loads/BC and create list of loads datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  listLoadPrep      : parameters to construct list of loads
! In  listLoad          : name of datastructure for list of loads
! In  jvBase            : JEVEUX base where to create objects
!
! --------------------------------------------------------------------------------------------------
!
    character(len=4), parameter :: phenom = "ACOU"
    integer(kind=8) :: iLoadList, iKeyword, indxLoadInList
    integer(kind=8) :: nbLoadList
    character(len=24) :: listLoad
    character(len=8), parameter :: funcCste = '&&ACDOCH'
    character(len=16) :: loadKeyword, loadCommand
    character(len=16), parameter :: loadApply = 'FIXE_CSTE'
    character(len=8) :: loadName, loadFunc
    character(len=13) :: loadPreObject
    aster_logical :: loadIsFunc, hasMultFunc
    character(len=8), pointer :: loadDble(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    listLoad = listLoadZ
    nbLoadList = 0
    indxLoadInList = 0
    iKeyword = 0

! - Get factor keyword to read
    call getLoadKeyword(loadKeyword)

! - Get number of loads for loads datastructure
    call getNbLoadsFromUser(phenom, loadKeyword, nbLoadList)

! - Create list of loads datastructure
    call creaListLoad(phenom, jvBase, nbLoadList, listLoad)

! - Create object to avoid same loads
    if (nbLoadList .ne. 0) then
        AS_ALLOCATE(vk8=loadDble, size=nbLoadList)
    end if

! - Add loads
    do iLoadList = 1, nbLoadList
! ----- Get access to load
        call getLoadIndexKeyword(loadKeyword, iKeyword)

! ----- Get current load
        call getLoadName(loadKeyword, iKeyword, &
                         iLoadList, nbLoadList, loadDble, &
                         loadName)

! ----- Get function applied to load
        call getLoadFunc(listLoadPrep, &
                         jvBase, funcCste, &
                         loadKeyword, iKeyword, &
                         loadFunc, hasMultFunc)

! ----- Get how to apply load
! ----- Only FIXE_CSTE

! ----- Get other load parameters
        call getLoadParameters(phenom, listLoadPrep%model, loadName, &
                               loadPreObject, loadCommand, loadIsFunc)
        ASSERT(.not. loadIsFunc)

! ----- Add acoustic loads
        call addLoadAcou(listLoad, &
                         loadName, loadFunc, &
                         loadApply, loadCommand, loadPreObject, &
                         loadIsFunc, &
                         indxLoadInList)
    end do

! - Debug
    !call listLoadDebug(listLoad)

! - Clean
    AS_DEALLOCATE(vk8=loadDble)
    call jedema()

end subroutine
