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
#include "asterf_types.h"
#include "contact_module.h"
!
interface
    subroutine laMatr(elem_dime   , l_axis        , nb_dofs, &
                    nb_lagr_c     , indi_lagc   , lagc_curr, &
                    gamma_c_nodes, &
                    nb_node_slav, elem_slav_code, slav_coor_init, slav_coor_curr,&
                    nb_node_mast, elem_mast_code, mast_coor_curr,&
                    proj_tole,  matr)
        integer, intent(in) :: elem_dime
        aster_logical, intent(in) :: l_axis
        integer, intent(in) :: nb_lagr_c, indi_lagc(9), nb_dofs
        character(len=8), intent(in) :: elem_slav_code, elem_mast_code
        integer, intent(in) :: nb_node_slav, nb_node_mast
        real(kind=8), intent(in) :: slav_coor_curr(3, 9), slav_coor_init(3,9)
        real(kind=8), intent(in) :: mast_coor_curr(3, 9)
        real(kind=8), intent(in) :: proj_tole, gamma_c_nodes(4), lagc_curr(4)
        real(kind=8), intent(inout) :: matr(MAX_CONT_DOFS, MAX_CONT_DOFS)
    end subroutine laMatr
end interface
