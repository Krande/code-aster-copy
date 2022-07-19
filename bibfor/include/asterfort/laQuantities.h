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
!
interface
    subroutine laQuantities(elem_dime, nb_node_slav, nb_node_mast, &
                        indi_lagc, &
                        slav_coor_init, mast_coor_init, &
                        slav_coor_curr, mast_coor_curr, &
                        slav_depl_curr, mast_depl_curr, lagc_curr)
        integer, intent(in) :: elem_dime, indi_lagc(9)
        integer, intent(in) :: nb_node_slav, nb_node_mast
        real(kind=8), intent(out) :: mast_coor_init(3, 9), slav_coor_init(3, 9)
        real(kind=8), intent(out) :: mast_coor_curr(3, 9), slav_coor_curr(3, 9)
        real(kind=8), intent(out) :: mast_depl_curr(3, 9), slav_depl_curr(3, 9)
        real(kind=8), intent(out) :: lagc_curr(4)
    end subroutine laQuantities
end interface
