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
    subroutine get_lagr_info(nume_equaz, i_equa, idx_gd, nb_node_lagr, list_node_lagr,&
                             nume_cmpz , ligrelz)
        character(len=*), intent(in) :: nume_equaz
        integer(kind=8), intent(in) :: i_equa
        integer(kind=8), intent(in) :: idx_gd
        integer(kind=8), intent(out) :: nb_node_lagr
        integer(kind=8), pointer :: list_node_lagr(:)
        integer(kind=8), optional, intent(out) :: nume_cmpz
        character(len=*), optional, intent(out) :: ligrelz
    end subroutine get_lagr_info
end interface
