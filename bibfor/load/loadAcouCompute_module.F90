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
! Module for the management of computation of loads (acoustic)
!
! ==================================================================================================
!
module loadAcouCompute_module
! ==================================================================================================
    use loadAcouCompute_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: isAcouLoadExist
    private :: getAcouNeumField
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/exisd.h"
#include "LoadTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getAcouNeumField
!
! Get name of input field to define acoustic Neumann loads
!
! In  indxNeuaType      : index of the type
! In  loadPreObject     : base JEVEUX name for object
! Out loadField         : name of input field to define load
!
! --------------------------------------------------------------------------------------------------
    subroutine getAcouNeumField(indxNeuaType, loadPreObjectZ, loadField)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeuaType
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=24), intent(out) :: loadField
!   ------------------------------------------------------------------------------------------------
!
        loadField = " "
        ASSERT(indxNeuaType .ge. 1 .and. indxNeuaType .le. LOAD_NEUA_NBTYPE)
        loadField = loadPreObjectZ(1:13)//acouLoadField(indxNeuaType)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isAcouLoadExist
!
! Detect existence of acoustic load
!
! In  indxNeuaType      : index of the type
! In  loadPreObject     : base JEVEUX name for object
! Out loadExist         : flag if load exists
! Out loadField         : standard input field
! In  hasInputField     : input field given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input field given by load (for EVOL_CHAR)
!
! --------------------------------------------------------------------------------------------------
    subroutine isAcouLoadExist(indxNeuaType, loadPreObjectZ, &
                               loadExist, &
                               loadField_, hasInputField_, inputLoadFieldZ_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeuaType
        character(len=*), intent(in) :: loadPreObjectZ
        aster_logical, intent(out) :: loadExist
        character(len=24), optional, intent(out) :: loadField_
        aster_logical, optional, intent(in) :: hasInputField_
        character(len=*), optional, intent(in) :: inputLoadFieldZ_
! ----- Local
        integer(kind=8) :: iret
        character(len=24) :: loadField, inputLoadField
        aster_logical :: hasInputField
!   ------------------------------------------------------------------------------------------------
!
        loadExist = ASTER_FALSE
        loadField = " "

! ----- Inputs
        hasInputField = ASTER_FALSE
        inputLoadField = " "
        if (present(hasInputField_)) then
            hasInputField = hasInputField_
        end if
        if (present(inputLoadFieldZ_)) then
            inputLoadField = inputLoadFieldZ_
        end if

! ----- Identify current load: get name of input field for this load
        if (hasInputField) then
            loadField = inputLoadField
        else
! --------- Field to detect
            if (indxNeuaType .le. LOAD_NEUT_NBTYPE) then
                call getAcouNeumField(indxNeuaType, loadPreObjectZ, loadField)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- Is this load exists ?
        iret = 0
        if (loadField .ne. " ") then
            call exisd('CHAMP_GD', loadField, iret)
        end if

! ----- Set outputs
        loadExist = iret .ne. 0
        if (present(loadField_)) then
            loadField_ = loadField
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module loadAcouCompute_module
