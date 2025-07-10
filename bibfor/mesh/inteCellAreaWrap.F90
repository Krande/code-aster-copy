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
subroutine inteCellAreaWrap(mesh, nbPoinInte, poinInteSlav, &
                            inteArea)
!
    use MeshPairing_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/mesh_pairing_type.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: nbPoinInte
    real(kind=8), intent(in) :: poinInteSlav(2, MAX_NB_INTE)
    real(kind=8), intent(out) :: inteArea
!
! --------------------------------------------------------------------------------------------------
!
! Pairing segment to segment
!
! Compute area of intersection
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : mesh
! In  nbPoinInte       : number of intersection points
! In  poinInte         : list of intersection points
! Out inteArea         : area of intersection
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: spaceDime
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('DIM_GEOM', mesh, 'MAILLAGE', repi=spaceDime)
    if (spaceDime .eq. 2) then
        spaceDime = 2
    else
        spaceDime = 3
    end if
    inteArea = 0.d0
    call inteCellArea(spaceDime, nbPoinInte, poinInteSlav, inteArea)
!
end subroutine
