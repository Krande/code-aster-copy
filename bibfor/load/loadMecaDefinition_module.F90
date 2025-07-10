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
! Module for the management of definition of loads (mechanics)
!
! ==================================================================================================
!
module loadMecaDefinition_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: creaLoadObje
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/wkvect.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! creaLoadObje
!
! Create general objects for definition of mechanical loads in AFFE_CHAR_MECA*
!
! In  jvBase            : base for JEVEUX
! In  load              : name of load
! In  modelZ            : name of model
!
! --------------------------------------------------------------------------------------------------
    subroutine creaLoadObje(jvBase, loadZ, modelZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(in) :: loadZ, modelZ
! ----- Locals
        character(len=8), pointer :: contNomo(:) => null(), contType(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        call wkvect(loadZ(1:8)//'.CHME.MODEL.NOMO', jvBase//' V K8', 1, vk8=contNomo)
        contNomo(1) = modelZ
        call wkvect(loadZ(1:8)//'.TYPE', jvBase//' V K8', 1, vk8=contType)
        contType(1) = 'MECA_RE'
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module loadMecaDefinition_module
