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
subroutine cmraff(mesh_in, mesh_out, level, info)
!
    use crea_maillage_module
!
    implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
!
    character(len=8), intent(in) :: mesh_in, mesh_out
    integer(kind=8), intent(in) :: level, info
!
! ----------------------------------------------------------------------
!         RAFFINEMENT UNIFORME DES MAILLES
! ----------------------------------------------------------------------
! IN        mesh_in   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  mesh_out  K8  NOM DU MAILLAGE TRANSFORME
! ----------------------------------------------------------------------
!
    type(Mmesh) :: mesh_raff
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Create new mesh
!
    call mesh_raff%init(mesh_in, info)
!
! - Refine
!
    call mesh_raff%refine(level)
!
! - Copy mesh
!
    call mesh_raff%copy_mesh(mesh_out)
!
! - Cleaning
!
    call mesh_raff%clean()
!
    call jedema()
end subroutine
