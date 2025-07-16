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
#include "asterf_types.h"

interface
    subroutine getelem(mesh, keywordfact, iocc, stop_void, list_elem, &
                       nb_elem, suffix, model, l_keep_propz, l_allz, onAllCells_)
        character(len=8), intent(in) :: mesh
        character(len=*), intent(in) :: keywordfact
        integer(kind=8), intent(in) :: iocc
        character(len=1), intent(in) :: stop_void
        integer(kind=8), intent(out) :: nb_elem
        character(len=24), intent(in) :: list_elem
        character(len=8), optional, intent(in) :: model
        character(len=*), optional, intent(in) :: suffix
        aster_logical, optional, intent(in) :: l_keep_propz, l_allz
        aster_logical, optional, intent(out) :: onAllCells_
    end subroutine getelem
end interface
