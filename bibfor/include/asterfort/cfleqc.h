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
    subroutine cfleqc(mesh       , sdcont_defi, nb_cont_zone, nb_cont_node, nb_cont_surf,&
                      v_poin_node, v_indi_node, nb_node_elim)
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: sdcont_defi
        integer(kind=8), intent(in) :: nb_cont_zone
        integer(kind=8), intent(in) :: nb_cont_surf
        integer(kind=8), intent(in) :: nb_cont_node
        integer(kind=8), pointer :: v_poin_node(:)
        integer(kind=8), pointer :: v_indi_node(:)
        integer(kind=8), intent(out) :: nb_node_elim
    end subroutine cfleqc
end interface
