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
#include "asterf_types.h"
!
interface
    subroutine aplcpgn(mesh, newgeo, &
                       mastConxInvName, mastNeighName, slavNeighName, &
                       pair_tole, dist_ratio, &
                       nb_elem_mast, list_elem_mast, nb_elem_slav, list_elem_slav, list_node_mast, &
                       nb_node_mast, meshPairing)
        use MeshPairing_module
        character(len=8), intent(in) :: mesh
        character(len=19), intent(in) :: newgeo
        character(len=24), intent(in) :: mastConxInvName
        character(len=24), intent(in) :: mastNeighName, slavNeighName
        real(kind=8), intent(in) :: pair_tole, dist_ratio
        integer(kind=8), intent(in) :: nb_elem_slav
        integer(kind=8), intent(in) :: nb_elem_mast
        integer(kind=8), intent(in) :: nb_node_mast
        integer(kind=8), intent(in) :: list_elem_mast(nb_elem_mast)
        integer(kind=8), intent(in) :: list_elem_slav(nb_elem_slav)
        integer(kind=8), intent(in) :: list_node_mast(nb_node_mast)
        type(MESH_PAIRING), intent(inout) :: meshPairing
    end subroutine aplcpgn
end interface
