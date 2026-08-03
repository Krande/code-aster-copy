! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
! Module for management of coordinates system
!
! ==================================================================================================
module coorSyst_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: setOrieFields, hasOrieField
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! setOrieFields
!
! Set orientation fields in input fields (computation)
!
! --------------------------------------------------------------------------------------------------
    subroutine setOrieFields(nbFieldInMax, lpain, lchin, &
                             nbFieldIn, caraElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbFieldInMax
        character(len=*), intent(inout) :: lpain(nbFieldInMax)
        character(len=*), intent(inout) :: lchin(nbFieldInMax)
        integer(kind=8), intent(inout) :: nbFieldIn
        character(len=*), optional, intent(in) :: caraElemZ
! ----- Local
        character(len=8) :: caraElem
        integer(kind=8) :: nbFieldAdd
!   ------------------------------------------------------------------------------------------------
!
        caraElem = caraElemZ
        nbFieldAdd = 3
        ASSERT(nbFieldIn+nbFieldAdd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PCAORIE'
        lchin(nbFieldIn+1) = caraElem(1:8)//'.CARORIEN'
        lpain(nbFieldIn+2) = 'PCACOQU'
        lchin(nbFieldIn+2) = caraElem(1:8)//'.CARCOQUE'
        lpain(nbFieldIn+3) = 'PCAMASS'
        lchin(nbFieldIn+3) = caraElem(1:8)//'.CARMASSI'
        nbFieldIn = nbFieldIn+nbFieldAdd
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! hasOrieField
!
! Detect orientation fields in input fields
!
! --------------------------------------------------------------------------------------------------
    function hasOrieField(jvCamass_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: hasOrieField
        integer(kind=8), optional, intent(out) :: jvCamass_
! ----- Local
        integer(kind=8) :: jvCamass
        integer(kind=8) :: iret, jtab(7)
!   ------------------------------------------------------------------------------------------------
!
        jvCamass = 0
        hasOrieField = ASTER_FALSE
        call tecach('NNO', 'PCAMASS', 'L', iret, nval=1, itab=jtab)
        if (iret .eq. 0) then
            jvCamass = jtab(1)
            hasOrieField = ASTER_TRUE
        else
            jvCamass = 0
            hasOrieField = ASTER_FALSE
        end if
        if (present(jvCamass_)) then
            jvCamass_ = jvCamass
        end if
!
!   ------------------------------------------------------------------------------------------------
    end function
!===================================================================================================
!===================================================================================================
!===================================================================================================
!===================================================================================================
!
end module coorSyst_module
