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
    subroutine laelem(nomte         , elem_dime     ,&
                    l_axis        , &
                    nb_dof        , nb_lagr_c       , indi_lagc   ,&
                    elem_slav_code, nb_node_slav,&
                    elem_mast_code, nb_node_mast)
        character(len=16), intent(in) :: nomte
        integer, intent(out) :: elem_dime
        aster_logical, intent(out) :: l_axis
        integer, intent(out) :: nb_dof
        integer, intent(out) :: nb_lagr_c
        integer, intent(out) :: indi_lagc(9)
        character(len=8), intent(out) :: elem_slav_code
        integer, intent(out) :: nb_node_slav
        character(len=8), intent(out) :: elem_mast_code
        integer, intent(out) :: nb_node_mast
    end subroutine laelem
end interface
