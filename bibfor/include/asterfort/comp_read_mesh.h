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
    subroutine comp_read_mesh(mesh          , keywordfact, iocc        ,&
                              list_elem_affe, l_affe_all , nb_elem_affe)
        character(len=8), intent(in) :: mesh
        character(len=16), intent(in) :: keywordfact
        integer(kind=8), intent(in) :: iocc
        character(len=24), intent(in) :: list_elem_affe
        aster_logical, intent(out) :: l_affe_all
        integer(kind=8), intent(out):: nb_elem_affe
    end subroutine comp_read_mesh
end interface
