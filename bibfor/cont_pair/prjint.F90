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
subroutine prjint(proj_tole, elem_dime, &
                  elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                  elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                  poin_inte, inte_weight, nb_poin_inte, &
                  inte_neigh_, ierror_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcodrm.h"
#include "asterfort/insema.h"
#include "asterfort/ptinma.h"
#include "asterfort/reereg.h"
#include "asterfort/apinte_norm.h"
#include "asterfort/apelem_getvertex.h"
#include "asterfort/apelem_inside.h"
#include "asterfort/apinte_weight.h"
#include "asterfort/apinte_chck.h"
#include "asterfort/apinte_prsl.h"
!
    real(kind=8), intent(in) :: proj_tole
    integer(kind=8), intent(in) :: elem_dime
    integer(kind=8), intent(in) :: elem_mast_nbnode
    real(kind=8), intent(in) :: elem_mast_coor(3, 9)
    character(len=8), intent(in) :: elem_mast_code
    integer(kind=8), intent(in) :: elem_slav_nbnode
    real(kind=8), intent(in) :: elem_slav_coor(3, 9)
    character(len=8), intent(in) :: elem_slav_code
    real(kind=8), intent(out) :: poin_inte(elem_dime-1, 16)
    real(kind=8), intent(out) :: inte_weight
    integer(kind=8), intent(out) :: nb_poin_inte
    integer(kind=8), optional, intent(inout) :: inte_neigh_(4)
    integer(kind=8), optional, intent(inout) :: ierror_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Pairing segment to segment
!
! Projection/intersection of elements in slave parametric space
!
! --------------------------------------------------------------------------------------------------
!
! In  proj_tole        : tolerance for projection
! In  elem_dime        : dimension of elements
! In  elem_mast_nbnode : number of nodes of master element
! In  elem_mast_coor   : coordinates of master element
! In  elem_mast_code   : code of master element
! In  elem_slav_nbnode : number of nodes for slave element
! In  elem_slav_coor   : coordinates of slave element
! In  elem_slav_code   : code of slave element
! Out poin_inte        : list (sorted) of intersection points
! Out inte_weight      : total weight of intersection
! Out nb_poin_inte     : number of intersection points
! IO  inte_neigh       : activation of neighbours of intersection
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: debug, l_inter
    real(kind=8) :: node_line_coop(elem_dime-1, 4)
    real(kind=8) :: proj_coop(elem_dime-1, 4), coor(3, 4)
    real(kind=8) :: xpt, ypt, xe(3), pt(3)
    real(kind=8) :: xp1, yp1, xp2, yp2
    integer(kind=8) :: test, list_next(16), nb_node_line
    integer(kind=8) :: i_node, iret
    character(len=8) :: elem_line_code
    real(kind=8) :: mast_norm(3), slav_norm(3)
    integer(kind=8) :: list_prev(16)
    integer(kind=8) :: inte_neigh(4)
!
! --------------------------------------------------------------------------------------------------
!
    nb_poin_inte = 0
    inte_weight = 0.d0
    iret = 0
    poin_inte(:, :) = 0.d0
    debug = ASTER_FALSE
    node_line_coop(elem_dime-1, 4) = 0.d0
    inte_neigh(1:4) = 0
    if (present(inte_neigh_)) then
        inte_neigh(:) = inte_neigh_(:)
    end if
    if (debug) then
        write (*, *) ". Projection/intersection"
    end if
!
! - Compute norms
!
    call apinte_norm(elem_dime, &
                     elem_mast_nbnode, elem_mast_coor, elem_mast_code, &
                     elem_slav_coor, elem_slav_code, &
                     mast_norm, slav_norm)
!
! - Project master nodes in slave element parametric space
!
    call apinte_prsl(proj_tole, elem_dime, &
                     elem_mast_nbnode, elem_mast_coor, &
                     elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                     proj_coop, iret)
    if (iret .eq. 1) then
        goto 99
    end if
!
! - Check if intersection is void or not
!
    call apinte_chck(proj_tole, elem_dime, &
                     elem_mast_nbnode, elem_mast_coor, &
                     elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                     proj_coop, mast_norm, slav_norm, &
                     l_inter)
    if (.not. l_inter) then
        goto 99
    end if
!
! - Get parametric coordinates of slave nodes (linear)
!
    call apelem_getvertex(elem_dime, elem_slav_code, &
                          node_line_coop, nb_node_line, elem_line_code, &
                          elem_slav_coor, proj_tole)
!
! - Set index of previous nodes
!
    do i_node = 2, nb_node_line
        list_prev(i_node) = i_node-1
    end do
    list_prev(1) = nb_node_line
!
! - Save projection of master nodes on slave element in list of intersection points
!
    call apelem_inside(proj_tole, elem_dime, elem_line_code, node_line_coop, &
                       elem_mast_nbnode, proj_coop, &
                       nb_poin_inte, poin_inte)
!
! - Add slave nodes if they are inside master element
!
    do i_node = 1, nb_node_line
! ----- Current coordinates of slave node
        xpt = node_line_coop(1, i_node)
        ypt = 0.d0
        if (elem_dime .eq. 3) then
            ypt = node_line_coop(2, i_node)
        end if
! ----- Test if point is inside element
        call ptinma(elem_mast_nbnode, elem_dime, elem_mast_code, proj_coop, proj_tole, &
                    xpt, ypt, test)
        if (test .eq. 1) then
            nb_poin_inte = nb_poin_inte+1
            poin_inte(1, nb_poin_inte) = xpt
            if (elem_dime .eq. 3) then
                poin_inte(2, nb_poin_inte) = ypt
            end if
            if (elem_dime .eq. 3) then
                inte_neigh(i_node) = 1
                inte_neigh(list_prev(i_node)) = 1
            else if (elem_dime .eq. 2) then
                inte_neigh(i_node) = 1
            else
                ASSERT(.false.)
            end if
        else if (test .eq. -1) then
            nb_poin_inte = 0
            poin_inte(1:elem_dime-1, 1:16) = 0.d0
            inte_neigh(1:4) = 0
            goto 99
        end if
    end do
!
! - Set index of next nodes
!
    do i_node = 2, elem_mast_nbnode
        list_next(i_node-1) = i_node
    end do
    list_next(elem_mast_nbnode) = 1
!
! - Intersection of edges
!
    if (elem_dime .eq. 3) then
        do i_node = 1, elem_mast_nbnode
!
! --------- Segment from edge of master element
!
            xp1 = proj_coop(1, i_node)
            yp1 = proj_coop(2, i_node)
            xp2 = proj_coop(1, list_next(i_node))
            yp2 = proj_coop(2, list_next(i_node))
!
! --------- Compute intersection between edge of master and slave element
!
            call insema(nb_node_line, elem_dime, node_line_coop, proj_tole, &
                        xp1, yp1, xp2, yp2, &
                        nb_poin_inte, poin_inte, inte_neigh)
        end do
        ASSERT(nb_poin_inte .le. 16)
    end if
!
! - Return in true parametric space for QUAD slave element
!
    if (elem_line_code .eq. "QU4") then
        do i_node = 1, 4
            coor(1:2, i_node) = node_line_coop(1:2, i_node)
            coor(3, i_node) = 0.d0
        end do
        do i_node = 1, nb_poin_inte
            pt(1) = poin_inte(1, i_node)
            pt(2) = poin_inte(2, i_node)
            pt(3) = 0.d0
            xe(:) = 0.d0
            call reereg('C', elem_line_code, 4, coor, &
                        pt, elem_dime, xe, test, proj_tole)
            if (test .eq. 1) then
                iret = 1
                nb_poin_inte = 0
                poin_inte(1:elem_dime-1, 1:16) = 0.d0
                inte_neigh(1:4) = 0
                goto 99
            end if
            poin_inte(1, i_node) = xe(1)
            poin_inte(2, i_node) = xe(2)
        end do
    end if

!
! - Sort list of intersection points
!
    if ((nb_poin_inte .gt. 2 .and. elem_dime .eq. 3) .or. &
        (nb_poin_inte .ge. 2 .and. elem_dime .eq. 2)) then
        call lcodrm(elem_dime, proj_tole, nb_poin_inte, poin_inte)
    end if
!
! - Check number of intersection point
!
    if (nb_poin_inte .gt. 8) then
        iret = 1
        nb_poin_inte = 0
        poin_inte(1:elem_dime-1, 1:16) = 0.d0
        inte_neigh(1:4) = 0
        goto 99
    end if
!
! - Compute weight of intersection
!
    call apinte_weight(elem_dime, nb_poin_inte, poin_inte, &
                       inte_weight)
!
! - No intersection exit
!
99  continue
    if (debug) then
        ASSERT(iret .eq. 0)
    end if
!
! - Copy
!
    if (present(inte_neigh_)) then
        inte_neigh_(1:4) = inte_neigh(1:4)
    end if
    if (present(ierror_)) then
        ierror_ = iret
    end if
!
end subroutine
