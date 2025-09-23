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

subroutine fix_mesh(mesh_in, mesh_out, info, remove_orphelan)
!
    use crea_maillage_module
!
    implicit none
#include "asterfort/cargeo.h"
#include "asterfort/infoma.h"
!
    character(len=8), intent(in) :: mesh_in, mesh_out
    integer(kind=8), intent(in) :: info, remove_orphelan
!
    type(Mmesh) :: mesh_conv
!
    call mesh_conv%init(mesh_in, info)
    call mesh_conv%fix(logical(remove_orphelan == 1, kind=1))
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
