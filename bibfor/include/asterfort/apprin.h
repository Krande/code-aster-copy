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
interface
    subroutine apprin(mesh          , newgeo        , pair_tole      ,nb_elem_mast  ,&
                      list_elem_mast, nb_elem_slav  , list_elem_slav ,elem_slav_flag,&
                      nb_mast_start , elem_mast_start,nb_slav_start  ,elem_slav_start)
        character(len=8), intent(in) :: mesh
        character(len=19), intent(in) :: newgeo
        real(kind=8), intent(in) :: pair_tole
        integer(kind=8), intent(in) :: nb_elem_mast
        integer(kind=8), intent(in) :: list_elem_mast(nb_elem_mast)
        integer(kind=8), intent(in) :: nb_elem_slav
        integer(kind=8), intent(in) :: list_elem_slav(nb_elem_slav)
        integer(kind=8), pointer :: elem_slav_flag(:)
        integer(kind=8), intent(out) :: nb_mast_start
        integer(kind=8), intent(out) :: elem_mast_start(nb_elem_slav)
        integer(kind=8), intent(out) :: nb_slav_start
        integer(kind=8), intent(out) :: elem_slav_start(nb_elem_slav)
    end subroutine apprin
end interface
