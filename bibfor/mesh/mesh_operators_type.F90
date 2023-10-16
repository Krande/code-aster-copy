! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
module mesh_operators_type
    !
    use Behaviour_type
    !
    implicit none
    !
#include "asterf_types.h"

! ==================================================================================================
!
! Global variables
!
! ==================================================================================================

! ==================================================================================================
!
! Derivated types - Parameters from MODI_MAILLAGE command
!
! ==================================================================================================

    type MESH_OPER_ORIE_SHELL
        aster_logical :: orieByVect = ASTER_FALSE
        real(kind=8) :: orieVect(3) = 0.d0

        integer :: nbGroupCell = 0
        character(len=24), pointer :: listOfGroupOfCell(:) => null()

        integer :: nodeNume = 0

    end type MESH_OPER_ORIE_SHELL

    type MESH_OPER_MODI_PARA
        integer :: orieShell = 0
        type(MESH_OPER_ORIE_SHELL), pointer :: meshOperOrieShell(:) => null()
        integer :: orieSkin = 0
        integer :: orieLine = 0
    end type MESH_OPER_MODI_PARA

!
end module mesh_operators_type
