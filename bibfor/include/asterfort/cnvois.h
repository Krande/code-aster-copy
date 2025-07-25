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
    subroutine cnvois(mesh      , list_elem , conx_inve , nb_elem, elem_indx_mini,&  
                    elem_indx_maxi, elem_neigh)
        character(len=8), intent(in) :: mesh
        character(len=24), intent(in) :: conx_inve
        integer(kind=8), intent(in) :: nb_elem
        integer(kind=8), intent(in) :: list_elem(nb_elem)
        integer(kind=8), intent(in) :: elem_indx_mini
        integer(kind=8), intent(in) :: elem_indx_maxi
        character(len=24), intent(in) :: elem_neigh
    end subroutine cnvois
end interface
