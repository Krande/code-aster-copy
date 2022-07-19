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
    subroutine laElemCont(elem_dime, coor_qp_sl, proj_tole, &
                    nb_node_slav, elem_slav_code, slav_coor_curr,&
                    nb_node_mast, elem_mast_code, mast_coor_curr,&
                    nb_lagr_c, lagc_curr, indi_lagc, gamma_c_nodes, hF, &
                    lagr_c, gap, gamma_c, projRmVal, l_cont_qp,&
                    dGap, d2Gap, mu_c)
        integer, intent(in) :: elem_dime
        integer, intent(in) :: nb_lagr_c, indi_lagc(9)
        character(len=8), intent(in) :: elem_slav_code, elem_mast_code
        integer, intent(in) :: nb_node_slav, nb_node_mast
        real(kind=8), intent(in) :: slav_coor_curr(3, 9), mast_coor_curr(3, 9)
        real(kind=8), intent(in) :: coor_qp_sl(2), hF
        real(kind=8), intent(in) :: proj_tole, gamma_c_nodes(4), lagc_curr(4)
        real(kind=8), intent(out) :: lagr_c, gap, gamma_c, projRmVal
        aster_logical, intent(out) :: l_cont_qp
        real(kind=8), intent(out), optional :: dGap(MAX_CONT_DOFS)
        real(kind=8), intent(out), optional :: d2Gap(MAX_CONT_DOFS, MAX_CONT_DOFS)
        real(kind=8), intent(out), optional :: mu_c(MAX_CONT_DOFS)
    end subroutine laElemCont
end interface

