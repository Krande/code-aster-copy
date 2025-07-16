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
! ==================================================================================================
!
! Module for the management of MGIS/MFRONT behaviours
!
! ==================================================================================================
!
module MGIS_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: hasMFront
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! hasMFront
!
! Detect MFront behaviour in compor maps
!
! --------------------------------------------------------------------------------------------------
    function hasMFront(comporMapZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: hasMFront
        character(len=*), intent(in) :: comporMapZ
! ----- Locals
        character(len=24) :: comporMap
        character(len=16), pointer :: comporVale(:) => null()
        integer(kind=8), pointer :: comporDesc(:) => null()
        integer(kind=8) :: comporValeSize, mapNbZone, iMapZone, mapNbCmpMax
        character(len=16) :: compor(COMPOR_SIZE), extern_addr
!   ------------------------------------------------------------------------------------------------
!
        hasMFront = ASTER_FALSE
        comporMap = comporMapZ

! ----- Acces to map
        call jeveuo(comporMap(1:19)//'.DESC', 'L', vi=comporDesc)
        call jeveuo(comporMap(1:19)//'.VALE', 'L', vk16=comporVale)
        call jelira(comporMap(1:19)//'.VALE', 'LONMAX', comporValeSize)
        mapNbZone = comporDesc(3)
        mapNbCmpMax = comporValeSize/comporDesc(2)

! ----- Detection
        do iMapZone = 1, mapNbZone
            compor = comporVale(iMapZone)
            extern_addr = comporVale(mapNbCmpMax*(iMapZone-1)+MGIS_ADDR)
            if (extern_addr .ne. "VIDE" .and. extern_addr .ne. " ") then
                hasMFront = ASTER_TRUE
                exit
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end function
!
end module MGIS_module
