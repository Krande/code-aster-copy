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
    subroutine gtvois(v_connex  , v_connex_lcum, list_elem, nb_elem   , elem_nume, elem_code,&
                      v_conx_inv, v_inv_lcum   , nb_neigh , list_neigh)
        integer(kind=8), pointer :: v_connex(:)
        integer(kind=8), pointer :: v_connex_lcum(:)
        integer(kind=8), pointer :: v_conx_inv(:)
        integer(kind=8), pointer :: v_inv_lcum(:)
        integer(kind=8), intent(in) :: nb_elem
        integer(kind=8), intent(in) :: list_elem(nb_elem)
        integer(kind=8), intent(in) :: elem_nume
        character(len=8), intent(in) :: elem_code
        integer(kind=8), intent(in) :: nb_neigh
        integer(kind=8), intent(out) :: list_neigh(4)
    end subroutine gtvois
end interface
