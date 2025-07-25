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
    subroutine drz13d(noma, ligrmo, type_vale, nb_node, list_node, &
                  cmp_index_dx, cmp_index_dy, cmp_index_dz, cmp_index_drx, cmp_index_dry,&
                  cmp_index_drz, lisrel, nom_noeuds)
        character(len=8), intent(in)  :: noma
        character(len=19), intent(in) :: ligrmo
        character(len=4), intent(in) :: type_vale
        integer(kind=8), intent(in) :: nb_node
        character(len=24), intent(in) :: list_node
        integer(kind=8), intent(in) :: cmp_index_dx
        integer(kind=8), intent(in) :: cmp_index_dy
        integer(kind=8), intent(in) :: cmp_index_dz
        integer(kind=8), intent(in) :: cmp_index_drx
        integer(kind=8), intent(in) :: cmp_index_dry
        integer(kind=8), intent(in) :: cmp_index_drz
        character(len=19), intent(in) :: lisrel
        character(len=8), intent(out) :: nom_noeuds(:)
    end subroutine drz13d
end interface
