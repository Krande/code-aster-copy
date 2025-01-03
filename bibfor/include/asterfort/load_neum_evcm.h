! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
    subroutine load_neum_evcm(inst_curr , load_name, i_load, ligrel_calc,&
                              nb_in_maxi, nb_in_prep, lpain    , lchin ,&
                              idx_matr , matr_elem)
        real(kind=8), intent(in) :: inst_curr
        character(len=8), intent(in) :: load_name
        integer, intent(in) :: i_load
        character(len=19), intent(in) :: ligrel_calc
        integer, intent(in) :: nb_in_maxi
        character(len=*), intent(inout) :: lpain(nb_in_maxi)
        character(len=*), intent(inout) :: lchin(nb_in_maxi)
        integer, intent(in) :: nb_in_prep
        integer, intent(inout) :: idx_matr
        character(len=19), intent(in) :: matr_elem
    end subroutine load_neum_evcm
end interface
