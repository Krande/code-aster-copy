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
!
#include "asterf_types.h"
!
interface
    subroutine getnode(mesh   , keywordfact, iocc, stop_void, list_node, &
                       nb_node, model      , suffix, elem_excl)
        character(len=8), intent(in) :: mesh
        character(len=16), intent(in) :: keywordfact
        integer(kind=8), intent(in) :: iocc
        character(len=1), intent(in) :: stop_void
        integer(kind=8), intent(out) :: nb_node
        character(len=24), intent(in) :: list_node
        character(len=8), intent(in), optional :: model
        character(len=*), intent(in), optional :: suffix
        aster_logical, intent(in), optional :: elem_excl
    end subroutine getnode
end interface
