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
subroutine intePoinCoorWrap(mesh, nodeCoorName, baseName, iPair, &
                            poinInteReal)
!
    use MeshPairing_module
    use mesh_type
    use mesh_cell_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mesh_pairing_type.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: nodeCoorName
    character(len=24), intent(in) :: baseName
    integer(kind=8), intent(in) :: iPair
    real(kind=8), intent(out)  :: poinInteReal(3, MAX_NB_INTE)
!
! --------------------------------------------------------------------------------------------------
!
! Pairing segment to segment
!
! Compute point of intersection in global space
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : mesh
! In  baseName         : JEVEUX base name for output objects
! In  nodeCoorName     : JEVEUX name for coordinates of nodes
! In  iPair            : index of pair
! Out poinInteReal     : coordinates of intersection points en global space
!
! --------------------------------------------------------------------------------------------------
!
    type(CELL_GEOM) :: cellSlav
    integer(kind=8) :: nbPoinInte
    real(kind=8) :: poinInteSlav(2, MAX_NB_INTE)
    real(kind=8), pointer :: nodeCoor(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    poinInteReal = 0.d0

! - Access to updated geometry
    call jeveuo(nodeCoorName(1:19)//'.VALE', 'L', vr=nodeCoor)

! - Access to pair objects
    call getPairJV(mesh, baseName, nodeCoor, iPair+1, cellSlav_=cellSlav)

! - Get coordinates of intersection points in slave parametric space
    call getInteJV(baseName, iPair+1, nbPoinInte, poinInteSlav)

! - Project coordinates
    call intePoinCoor(cellSlav, nbPoinInte, poinInteSlav, poinInteReal)
!
end subroutine
