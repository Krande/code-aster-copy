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
interface
    subroutine pj2dco(typeSelect, &
                      entity1, entity2, &
                      nbCellSelect1, listCellSelect1, &
                      nbNodeSelect2, listNodeSelect2, &
                      geom1, geom2, corrMesh, &
                      l_dmax, dmax, dala, &
                      listInterc_, nbInterc_)
        character(len=*), intent(in) :: typeSelect
        character(len=8), intent(in) :: entity1, entity2
        integer(kind=8), intent(in) :: nbCellSelect1, listCellSelect1(*)
        integer(kind=8), intent(in) :: nbNodeSelect2, listNodeSelect2(*)
        character(len=*), intent(in) :: geom1, geom2
        character(len=16), intent(in)  :: corrMesh
        aster_logical, intent(in) :: l_dmax
        real(kind=8), intent(in) :: dmax, dala
        character(len=16), optional, intent(in)  :: listInterc_
        integer(kind=8), optional, intent(in)  :: nbInterc_
    end subroutine pj2dco
end interface
