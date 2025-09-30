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

subroutine fix_mesh(mesh_in, mesh_out, remove_orphan, positive_volume, outward_normal, &
                    double_nodes, double_cells, tole, info)
!
    use crea_maillage_module
!
    implicit none
#include "asterf_types.h"
#include "asterfort/cargeo.h"
#include "asterfort/infoma.h"
!
    character(len=8), intent(in) :: mesh_in, mesh_out
    integer(kind=8), intent(in) :: info, remove_orphan, positive_volume, outward_normal
    integer(kind=8), intent(in) :: double_nodes, double_cells
    real(kind=8), intent(in) :: tole
!
    type(Mmesh) :: mesh_conv
    aster_logical :: ro, on, pv, dn, dc
!
    ro = remove_orphan == 1
    on = outward_normal == 1
    pv = positive_volume == 1
    dn = double_nodes == 1
    dc = double_cells == 1
!
    call mesh_conv%init(mesh_in, info, convert_max=ASTER_FALSE)
    call mesh_conv%fix(ro, on, pv, dn, dc, tole)
    call mesh_conv%check_mesh()
    call mesh_conv%copy_mesh(mesh_out)
    call mesh_conv%clean()
!
! - Update parameters for modified mesh (bounding box and dimensions)
    call cargeo(mesh_out)

! - Verbose
    call infoma(mesh_out, info)
!
end subroutine
