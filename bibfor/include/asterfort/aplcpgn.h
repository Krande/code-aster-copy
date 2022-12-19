! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
    subroutine aplcpgn(mesh, newgeo, zone,  pair_method,  pair_tole, dist_ratio, &
                      nb_elem_mast, list_elem_mast, nb_elem_slav, list_elem_slav, list_node_mast,&
                      nb_node_mast , nb_pair_zone, list_pair_zone, list_nbptit_zone,&
                      list_ptitsl_zone)
        character(len=8), intent(in) :: mesh
        character(len=19), intent(in) :: newgeo
        character(len=19), intent(in) :: zone
        real(kind=8), intent(in) :: pair_tole, dist_ratio
        integer, intent(in) :: nb_elem_slav
        integer, intent(in) :: nb_elem_mast
        integer, intent(in) :: nb_node_mast
        integer, intent(in) :: list_elem_mast(nb_elem_mast)
        integer, intent(in) :: list_elem_slav(nb_elem_slav)
        integer, intent(in) :: list_node_mast(nb_node_mast)
        integer, intent(out) :: nb_pair_zone
        character(len=19), intent(in) :: list_pair_zone, list_nbptit_zone
        character(len=19), intent(in) :: list_ptitsl_zone
        character(len=24), intent(in) :: pair_method
    end subroutine aplcpgn
end interface
