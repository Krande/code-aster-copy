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
!
subroutine interNodesEdge(proj_tole       , elem_dime     , &
                       elem_mast_code, elem_slave_code, &
                       proj_coor       , nb_node_proj, &
                       nb_poin_inte,  poin_inte, inte_neigh)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/apelem_getvertex_n.h"
#include "asterfort/insema.h"

!
real(kind=8), intent(in) :: proj_tole
integer, intent(in) :: elem_dime
character(len=8), intent(in) :: elem_mast_code, elem_slave_code
real(kind=8), intent(in) :: proj_coor(elem_dime-1,9)
integer, intent(in) :: nb_node_proj
integer, intent(inout) :: inte_neigh(4), nb_poin_inte
real(kind=8), intent(inout) :: poin_inte(elem_dime-1,16)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Intersection of edges in 3D
!
! --------------------------------------------------------------------------------------------------
!
! In  proj_tole        : tolerance for projection
! In  elem_dime        : dimension of elements
! In  elem_mast_nbnode : number of nodes of master element
! In  elem_slav_nbnode : number of nodes for slave element
! In  elem_slav_coor   : coordinates of slave element
! In  elem_mast_code   : code of master element
! Out proj_coor        : projection of slave nodes on master element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: elem_mast_line_coop(elem_dime-1,4)
    integer :: elem_mast_line_nbnode, i_node
    integer :: list_next(8)
    character(len=8) :: elem_mast_line_code
    real(kind=8) :: xp1, yp1, xp2, yp2
!
! - Compute intersection in parametric master space
!
    if(elem_dime == 3) then
        ASSERT(nb_node_proj <= 8)
!
! - Get parametric coordinates of master nodes (linear)
!
        call apelem_getvertex_n(elem_dime, elem_mast_code,&
                                elem_mast_line_coop, elem_mast_line_nbnode, &
                                elem_mast_line_code)
!
! - Set index of next nodes
!
        list_next = 0
        if(elem_slave_code == "TR3") then
            list_next(1:3) = [2, 3, 1]
        elseif(elem_slave_code == "TR6") then
            list_next(1:6) = [4, 2, 5, 3, 6, 1]
        elseif(elem_slave_code == "QU4") then
            list_next(1:4) = [2, 3, 4, 1]
        elseif(elem_slave_code == "QU8" .or. elem_slave_code == "QU9") then
            list_next(1:8) = [5, 2, 6, 3, 7, 4, 8, 1]
        else
            ASSERT(ASTER_FALSE)
        end if
!
! - Intersection of edges
!
        do i_node = 1, nb_node_proj
!
! --------- Segment from edge of projected slave cell
!
            xp1 = proj_coor(1,i_node)
            yp1 = proj_coor(2,i_node)
            xp2 = proj_coor(1,list_next(i_node))
            yp2 = proj_coor(2,list_next(i_node))
!
! --------- Compute intersection between edge of master and projected slave cells
!
            call insema(elem_mast_line_nbnode, elem_dime, elem_mast_line_coop, proj_tole,&
                        xp1, yp1, xp2, yp2,&
                        nb_poin_inte, poin_inte, inte_neigh)
        end do
        inte_neigh = 1
    end if
!
end subroutine
