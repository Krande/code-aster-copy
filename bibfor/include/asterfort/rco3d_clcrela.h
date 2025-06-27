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
    subroutine rco3d_clcrela(ligrel, noma, nb_pairs, nbnocot,&
        list_total_no_co, map_noco_pair, map_noco_nbelem, &
        map_noco_nbnoco, resuelem, fonrez, lisrel )
        character(len=19), intent(in) :: ligrel, resuelem, lisrel
        character(len=8), intent(in) :: noma
        integer(kind=8), intent(in) :: nb_pairs, nbnocot
        integer(kind=8), intent(in) :: map_noco_pair(:,:,:)
        integer(kind=8), intent(in) :: map_noco_nbnoco(:,:,:)
        integer(kind=8), intent(in) :: map_noco_nbelem(:,:)
        integer(kind=8), pointer, intent(in) :: list_total_no_co(:)
        character(len=*), intent(in) :: fonrez
    end subroutine rco3d_clcrela
end interface
