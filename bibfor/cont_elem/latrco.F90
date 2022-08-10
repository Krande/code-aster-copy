! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine latrco(i_tria, nb_poin_inte, poin_inte, tria_coor)
!
implicit none
!
#include "asterfort/assert.h"
!
!
    integer, intent(in) :: i_tria
    integer, intent(in) :: nb_poin_inte
    real(kind=8), intent(in) :: poin_inte(2, 8)
    real(kind=8), intent(out) :: tria_coor(2, 3)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (LAC) - Elementary computations
!
! Coordinates of current triangle
!
! --------------------------------------------------------------------------------------------------
!
! In  i_tria           : index of current triangle
! In  tria_node        : list of triangles (defined by index of intersection points)
! In  poin_inte        : list (sorted) of intersection points
! Out tria_coor        : coordinates of current triangle
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i_node, i_node1, i_node2
    real(kind=8) ::  barycenter(2)
!
! --------------------------------------------------------------------------------------------------
!
    if(nb_poin_inte == 3) then
        tria_coor(1:2, 1:3) = poin_inte(1:2,1:3)
    else
        barycenter = 0.d0
        do i_node=1, nb_poin_inte
            barycenter(1:2) = barycenter(1:2) + poin_inte(1:2,i_node)
        end do
        barycenter = barycenter / real(nb_poin_inte, kind=8)
!
        i_node1 = i_tria
        i_node2 = i_tria + 1
        if(i_tria == nb_poin_inte) i_node2 = 1

        tria_coor(1:2,1) = poin_inte(1:2,i_node1)
        tria_coor(1:2,2) = poin_inte(1:2,i_node2)
        tria_coor(1:2,3) = barycenter
    end if
!
end subroutine
