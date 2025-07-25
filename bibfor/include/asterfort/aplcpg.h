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
!
#include "asterf_types.h"
!
interface
    subroutine aplcpg(mesh        , newgeo        , sdappa      , i_zone        , pair_tole,&
                      nb_elem_mast, list_elem_mast, nb_elem_slav, list_elem_slav, &
                      nb_pair_zone, list_pair_zone, list_nbptit_zone, list_ptitsl_zone,&
                      list_ptitma_zone,li_ptgausma_zone,i_proc      , nb_proc, pair_method)
        character(len=8), intent(in) :: mesh
        character(len=19), intent(in) :: newgeo
        character(len=19), intent(in) :: sdappa
        integer(kind=8), intent(in) :: i_zone
        real(kind=8), intent(in) :: pair_tole
        integer(kind=8), intent(in) :: nb_elem_slav
        integer(kind=8), intent(in) :: nb_elem_mast
        integer(kind=8), intent(in) :: list_elem_mast(nb_elem_mast)
        integer(kind=8), intent(in) :: list_elem_slav(nb_elem_slav)
        integer(kind=8), intent(inout) :: nb_pair_zone
        integer(kind=8), pointer :: list_pair_zone(:)
        integer(kind=8), pointer :: list_nbptit_zone(:)
        real(kind=8), pointer :: list_ptitsl_zone(:)
        real(kind=8), pointer :: list_ptitma_zone(:)
        real(kind=8), pointer :: li_ptgausma_zone(:)
        integer(kind=8), intent(in) :: i_proc
        integer(kind=8), intent(in) :: nb_proc
        character(len=24), intent(in) :: pair_method
    end subroutine aplcpg
end interface
