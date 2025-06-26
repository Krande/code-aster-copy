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
subroutine minMaxEdgesWrap(meshZ, groupNameZ, edgeMin, edgeMax)
!
    use mesh_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: meshZ, groupNameZ
    real(kind=8), intent(out) :: edgeMin, edgeMax
!
! --------------------------------------------------------------------------------------------------
!
! Mesh
!
! Compute min and max edges on group of cells
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : mesh
! In  groupNameZ       : name of cell group
! Out edgeMin          : length of smallest edge
! Out edgeMax          : length of greatest edge
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbCell
    integer(kind=8), pointer :: listCellNume(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call getCellsFromGroup(meshZ, groupNameZ, nbCell, listCellNume)
    if (nbCell .eq. 0) then
        call utmess('F', 'MESH3_13', sk=groupNameZ(1:24))
    else
        call compMinMaxEdges(meshZ, edgeMin, edgeMax, &
                             nbCell, listCellNume)
    end if
!
end subroutine
