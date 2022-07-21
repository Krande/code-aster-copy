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
    subroutine mmelem_data_laga(l_axi_         , &
                             typg_slav_name_, typg_mast_name_,&
                             nb_cont_type_  , nb_node_elem_  ,&
                             typg_cont_nume_, &
                             typf_cont_nume_, &
                             typf_frot_nume_, &
                             get_elem_indx_)
        aster_logical, intent(in) :: l_axi_
        character(len=8), intent(in) :: typg_slav_name_
        character(len=8), intent(in) :: typg_mast_name_
        integer, intent(out) :: nb_cont_type_
        integer, intent(out) :: nb_node_elem_
        integer, intent(out) :: typg_cont_nume_
        integer, intent(out) :: typf_cont_nume_
        integer, intent(out) :: typf_frot_nume_
        integer, intent(out) :: get_elem_indx_
    end subroutine mmelem_data_laga
end interface
