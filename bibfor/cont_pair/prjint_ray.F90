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
    subroutine prjint_ray(proj_tole       , elem_dime     , &
                  elem_mast_nbnode, elem_mast_coor, elem_mast_code,&
                  elem_slav_nbnode, elem_slav_coor, elem_slav_code,&
                  poin_inte       , inte_weight   , nb_poin_inte  ,&
                  inte_neigh_, ierror_)
!
implicit none
!
#include "asterf_types.h"


#include "asterfort/apinte_prsl.h"
#include "asterfort/apinte_weight.h"
#include "asterfort/assert.h"
#include "asterfort/cfadju.h"
#include "asterfort/lcodrm.h"
#include "asterfort/reereg.h"
#include "asterfort/reerel.h"
#include "asterfort/projMaAndCheck.h"
#include "asterfort/interNodesInside.h"
#include "asterfort/interNodesEdge.h"
#include "asterfort/dctest.h"
!
real(kind=8), intent(in) :: proj_tole
integer, intent(in) :: elem_dime
integer, intent(in) :: elem_mast_nbnode
real(kind=8), intent(in) :: elem_mast_coor(3,9)
character(len=8), intent(in) :: elem_mast_code
integer, intent(in) :: elem_slav_nbnode
real(kind=8), intent(in) :: elem_slav_coor(3,9)
character(len=8), intent(in) :: elem_slav_code
real(kind=8), intent(out) :: poin_inte(elem_dime-1,8)
real(kind=8), intent(out) :: inte_weight
integer, intent(out) :: nb_poin_inte
integer, optional, intent(inout) :: inte_neigh_(4)
integer, optional, intent(inout) :: ierror_
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
! In  elem_mast_nbnode : number of nodes of master cell
! In  elem_mast_coor   : coordinates of master cell
! In  elem_mast_code   : code of master cell
! In  elem_slav_nbnode : number of nodes for slave cell
! In  elem_slav_coor   : coordinates of slave cell
! In  elem_slav_code   : code of slave cell
! Out poin_inte        : list (sorted) of intersection points
! Out inte_weight      : total weight of intersection
! Out nb_poin_inte     : number of intersection points
! IO  inte_neigh       : activation of neighbours of intersection
! Out iret : 0 - Ok and intersection is not void
!            1 - Ok and intersection is void
!            2 - Error during projection
! --------------------------------------------------------------------------------------------------
!
    aster_logical, parameter :: debug = ASTER_FALSE
    aster_logical :: error
    real(kind=8) :: proj_coop(elem_dime-1,9), coor_inte_ma(3, 9)
    real(kind=8) :: inte_weight_ma, xe(3)
    integer :: i_node, iret, nb_node_proj, elem_slav_line_nbnode
    character(len=8) :: elem_slav_line_code, elem_mast_line_code
    integer :: elem_mast_line_nbnode, elin_nbnode(1), test
    integer :: inte_neigh(4), elin_sub(1,4), elin_nbsub, nb_poin_inte_ma
    real(kind=8) :: poin_inte_ma(elem_dime-1,16)
!
! --------------------------------------------------------------------------------------------------
!
    nb_poin_inte       = 0
    inte_weight        = 0.d0
    iret               = 0
    poin_inte          = 0.d0
    poin_inte_ma       = 0.d0
    inte_neigh = 0
    error = ASTER_FALSE
!
    if (debug) then
        write(*,*) ". Projection/intersection with raytracing"
    endif
    ! print*, "MA: ", elem_mast_coor(1:3,1:elem_mast_nbnode)
    ! print*, "ES: ", elem_slav_coor(1:3,1:elem_slav_nbnode)
!
! - Linearized master and slave cell
!
    call dctest(elem_mast_code, elin_sub, elin_nbnode, elin_nbsub, elem_mast_line_code)
    ASSERT(elin_nbsub == 1)
    elem_mast_line_nbnode = elin_nbnode(1)
    call dctest(elem_slav_code, elin_sub, elin_nbnode, elin_nbsub, elem_slav_line_code)
    ASSERT(elin_nbsub == 1)
    elem_slav_line_nbnode = elin_nbnode(1)
!
! - Quick test on linearized master cell to detect empty intersection
!   Projection on linear master cell
!
    call projMaAndCheck(proj_tole, elem_dime, &
                        elem_mast_line_nbnode, elem_mast_coor, elem_mast_line_code,&
                        elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                        proj_coop, nb_node_proj, iret)
!
    ASSERT(iret < 2)
    if(iret == 1) then
        if (debug) then
            write(*,*) ".. Projection/intersection is empty (quick check 1)"
        endif
        go to 99
    end if
!
    !call interNodesInside(proj_tole, elem_dime, &
    !                      elem_mast_line_code, elem_slav_code, &
    !                      proj_coop, nb_node_proj, nb_poin_inte_ma, poin_inte_ma, inte_neigh)
!
! - If intersection is empty -> exit
!
    !if(nb_poin_inte_ma == 0) then
    !    if (debug) then
    !        write(*,*) ".. Projection/intersection is empty (quick check)"
    !    endif
    !    go to 99
    !end if
!
! - Intersection is not empty
! - Projection on master cell
!
    call projMaAndCheck(proj_tole, elem_dime, &
                        elem_mast_nbnode, elem_mast_coor, elem_mast_code,&
                        elem_slav_nbnode, elem_slav_coor, elem_slav_code, &
                        proj_coop, nb_node_proj, iret)
!
    if(iret == 2) then
        error = ASTER_TRUE
        goto 100
    end if
!
! - Compute intersection in parametric master space
!
    call interNodesInside(proj_tole, elem_dime, &
                        elem_mast_code, elem_slav_code, &
                        proj_coop, nb_node_proj, nb_poin_inte_ma, poin_inte_ma, inte_neigh)
!
    call interNodesEdge(proj_tole, elem_dime, &
                        elem_mast_code, elem_slav_code, &
                        proj_coop, nb_node_proj, nb_poin_inte_ma, poin_inte_ma, inte_neigh)
!
    if(nb_poin_inte_ma == 0) then
        error = ASTER_TRUE
        goto 100
    end if
!
    ASSERT(nb_poin_inte_ma .le. 16)
!
! - Sort list of intersection points
!
    if ((nb_poin_inte_ma .gt. 2 .and. elem_dime == 3) .or.&
        (nb_poin_inte_ma .ge. 2 .and. elem_dime == 2)) then
        call lcodrm(elem_dime, proj_tole, nb_poin_inte_ma, poin_inte_ma)
    endif
!
    if(nb_poin_inte_ma == 0 .or. nb_poin_inte_ma > 8) then
        iret = 1
        go to 99
    end if
!
! - Compute weight of intersection
!
    call apinte_weight(elem_dime, nb_poin_inte_ma, poin_inte_ma, inte_weight_ma)
!
    if(inte_weight_ma <= proj_tole) then
        iret = 1
        go to 99
    end if
!
! - Return in real master space
!
    coor_inte_ma = 0
    do i_node = 1, nb_poin_inte_ma
        xe = 0.d0
        xe(1:2) = poin_inte_ma(1:2, i_node)
        call reerel(elem_mast_code, elem_mast_nbnode, 3, elem_mast_coor, xe, &
                    coor_inte_ma(1:3, i_node))
    end do
!
! - Return in parametric slave space using orthogonal projection
!
    call apinte_prsl(proj_tole    , elem_dime     , &
                     nb_poin_inte_ma , coor_inte_ma, &
                     elem_slav_nbnode, elem_slav_coor, elem_slav_code,&
                     poin_inte       , iret)
    if(iret == 1) then
        error = ASTER_TRUE
        go to 100
    end if
!
! - All nodes have to be inside slace cell
!
    nb_poin_inte = nb_poin_inte_ma
    do i_node = 1, nb_poin_inte
! ----- Test if point is inside element
        call cfadju(elem_slav_line_code, poin_inte(1,i_node), poin_inte(2,i_node), &
                    proj_tole, test)

        if(test == 2) then
            error = ASTER_TRUE
            go to 100
        end if
    end do
!
! - Sort list of intersection points
!
    if ((nb_poin_inte .gt. 2 .and. elem_dime == 3) .or.&
        (nb_poin_inte .ge. 2 .and. elem_dime == 2)) then
        call lcodrm(elem_dime, proj_tole, nb_poin_inte, poin_inte)
    endif
    ASSERT(nb_poin_inte <= 8)
!
! - Compute weight of intersection
!
    call apinte_weight(elem_dime, nb_poin_inte, poin_inte, inte_weight)
!
! - Error
!
    100 continue
    if (error) then
        iret = 2
        nb_poin_inte = 0
        poin_inte = 0.d0
        inte_neigh = 0
        inte_weight = 0.d0
        !ASSERT(ASTER_FALSE)
    end if
!
! - No intersection exit
!
99  continue
    if (debug) then
        ASSERT(iret < 2)
    end if
!
! - Copy
!
    if (present(inte_neigh_)) then
        inte_neigh_(1:4) = inte_neigh(1:4)
    endif
    if (present(ierror_)) then
        ierror_=iret
    endif
!
end subroutine
