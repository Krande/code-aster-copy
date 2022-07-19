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
    subroutine getQuadCont(elem_dime, l_axis, &
                        nb_node_slav, elem_slav_code, elem_slav_coor,&
                        elem_mast_code, &
                        nb_qp, coor_qp     , weight_qp )
        integer, intent(in) :: elem_dime, nb_node_slav
        aster_logical, intent(in) :: l_axis
        character(len=8), intent(in) :: elem_slav_code, elem_mast_code
        real(kind=8), intent(in) :: elem_slav_coor(3, 9)
        real(kind=8), intent(out) :: coor_qp(2, 48), weight_qp(48)
        integer, intent(out) :: nb_qp
    end subroutine getQuadCont
end interface
