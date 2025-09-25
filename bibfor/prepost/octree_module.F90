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
module octree_module
!
    implicit none
!
    private
!
#include "asterf_types.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
!
    integer(kind=8), parameter, private :: nb_div = 3
!
    type :: OctreeNode
        real(kind=8) :: xmin(3) = 0.d0, dx(3) = 1.0
        type(OctreeNode), pointer :: children(:, :, :) => null()
        integer(kind=8), allocatable :: list_pts(:)
        real(kind=8), allocatable :: coordinates(:, :)
        integer(kind=8) :: nb_pts = 0
! ----- member functions
    contains
        procedure, public, pass :: split => split_node
        procedure, public, pass :: free => free_node
        procedure, public, pass :: is_strictly_inside
        procedure, public, pass :: find_children
    end type OctreeNode
!
    type Octree
        real(kind=8) :: xmin(3) = 0.d0, xmax(3) = 0.d0
        integer(kind=8) :: max_pts_by_node = 0, dim = 3, nb_pt = 0
        integer(kind=8), allocatable :: list_pts(:)
        type(OctreeNode) :: node
! ----- member functions
    contains
        procedure, public, pass :: init => init_octree
        procedure, public, pass :: free => free_octree
        procedure, public, pass :: get_pts_around
    end type
!
!===================================================================================================
!
!===================================================================================================
!
    public :: Octree, OctreeNode
    private :: split_node, free_node, init_octree, free_octree, get_pts_around
!
contains
!
    recursive subroutine free_node(this)
        class(OctreeNode), intent(inout) :: this
!
        integer(kind=8) :: i, j, k
!
        if (allocated(this%list_pts)) then
            this%nb_pts = 0
            deallocate (this%list_pts)
            deallocate (this%coordinates)
        else
            if (associated(this%children)) then
                do i = 1, size(this%children, 1)
                    do j = 1, size(this%children, 2)
                        do k = 1, size(this%children, 3)
                            call this%children(i, j, k)%free()
                        end do
                    end do
                end do
                deallocate (this%children)
            end if
        end if
!
    end subroutine free_node
!
! ==================================================================================================
!
!
    function is_strictly_inside(this, coor, distance) result(inside)
        class(OctreeNode), intent(in) :: this
        real(kind=8), intent(in) :: distance, coor(3)
!
        aster_logical :: inside
        inside = ASTER_FALSE
!
        if (this%xmin(1) <= (coor(1)-distance) .and. &
            (coor(1)+distance) <= (this%xmin(1)+this%dx(1))) then
            if (this%xmin(2) <= (coor(2)-distance) .and. &
                (coor(2)+distance) <= (this%xmin(2)+this%dx(2))) then
                if (size(this%children, 3) == 3) then
                    if (this%xmin(3) <= (coor(3)-distance) .and. &
                        (coor(3)+distance) <= (this%xmin(3)+this%dx(3))) then
                        inside = ASTER_TRUE
                    end if
                else
                    inside = ASTER_TRUE
                end if
            end if
        end if
!
    end function
!
! ==================================================================================================
!
    subroutine find_children(this, coor, child_id)
        class(OctreeNode), intent(in) :: this
        real(kind=8), intent(in) :: coor(3)
        integer(kind=8), intent(inout) :: child_id(3)
!
        integer(kind=8) :: i_dim, dim
!
        dim = 3
        if (size(this%children, 3) < nb_div) then
            dim = 2
        end if
        child_id = 1
        do i_dim = 1, dim
            child_id(i_dim) = min(nb_div, max(1, ceiling((coor(i_dim)-this%xmin(i_dim)) &
                                                         /this%dx(i_dim), kind=8)))
        end do
!
    end subroutine
!
! ==================================================================================================
!
    subroutine free_octree(this)
        class(Octree), intent(inout) :: this
!
        call this%node%free()
!
    end subroutine free_octree
!
! ==================================================================================================
!
    subroutine init_octree(this, nb_pts, coordinates, max_nb_pts)
        class(Octree), intent(inout) :: this
        integer(kind=8), intent(in) :: nb_pts, max_nb_pts
        real(kind=8), intent(in) :: coordinates(3, nb_pts)
!
        integer(kind=8), allocatable :: sub_pts(:)
        integer(kind=8) :: i_pt, i_dim
!
        this%max_pts_by_node = max_nb_pts
!
        this%xmin = R8MAEM()
        this%xmax = -R8MAEM()
!
        do i_pt = 1, nb_pts
            do i_dim = 1, 3
                this%xmin(i_dim) = min(this%xmin(i_dim), coordinates(i_dim, i_pt))
                this%xmax(i_dim) = max(this%xmax(i_dim), coordinates(i_dim, i_pt))
            end do
        end do
!
        this%dim = 3
        if (abs(this%xmax(3)-this%xmin(3)) < 1d-12) then
            this%dim = 2
        end if
!
        this%node%xmin = this%xmin
        this%node%dx = 1.0d0
        do i_dim = 1, this%dim
            this%node%dx = (this%xmax(i_dim)-this%node%xmin(i_dim))/real(nb_div, kind=8)
        end do
!
        allocate (sub_pts(nb_pts))
!
        do i_pt = 1, nb_pts
            sub_pts(i_pt) = i_pt
        end do
!
        call this%node%split(nb_pts, coordinates, this%dim, max_nb_pts, nb_pts, sub_pts)
        deallocate (sub_pts)
!
    end subroutine init_octree
!
! ==================================================================================================
!
    subroutine split_node(this, nb_pts, coordinates, dim, max_nb_pts, &
                          nb_sub_pt, list_sub_pt)
        class(OctreeNode), intent(inout) :: this
        integer(kind=8), intent(in) :: nb_pts, dim, nb_sub_pt, list_sub_pt(nb_sub_pt), max_nb_pts
        real(kind=8), intent(in) :: coordinates(3, nb_pts)
!
        integer(kind=8), allocatable :: sub_pts(:, :, :, :)
        integer(kind=8) :: i_pt, oi(3), pt_id, i, j, k, nb_sub_pts(3, 3, 3), kmax
!
        this%nb_pts = nb_sub_pt
        if (nb_sub_pt <= max_nb_pts) then
            if (nb_sub_pt > 0) then
                allocate (this%list_pts(nb_sub_pt))
                allocate (this%coordinates(3, nb_sub_pt))
!
                this%list_pts(1:nb_sub_pt) = list_sub_pt(1:nb_sub_pt)
                do i_pt = 1, nb_sub_pt
                    this%coordinates(1:3, i_pt) = coordinates(1:3, list_sub_pt(i_pt))
                end do
            end if
        else
!
            kmax = nb_div
            if (dim < 3) then
                kmax = 1
            end if
!
            allocate (this%children(nb_div, nb_div, kmax))
!
            do i = 1, nb_div
                do j = 1, nb_div
                    do k = 1, kmax
                        this%children(i, j, k)%dx(1) = this%dx(1)/real(nb_div, kind=8)
                        this%children(i, j, k)%dx(2) = this%dx(2)/real(nb_div, kind=8)
                        this%children(i, j, k)%dx(3) = this%dx(3)/real(kmax, kind=8)
                        this%children(i, j, k)%xmin(1) = this%xmin(1)+ &
                                                         (i-1)*this%children(i, j, k)%dx(1)
                        this%children(i, j, k)%xmin(2) = this%xmin(2)+ &
                                                         (j-1)*this%children(i, j, k)%dx(2)
                        this%children(i, j, k)%xmin(3) = this%xmin(3)+ &
                                                         (k-1)*this%children(i, j, k)%dx(3)
                    end do
                end do
            end do
!
            allocate (sub_pts(nb_div, nb_div, kmax, nb_sub_pt))
            nb_sub_pts = 0
!
            do i_pt = 1, nb_sub_pt
                pt_id = list_sub_pt(i_pt)
                call this%find_children(coordinates(:, pt_id), oi)
                nb_sub_pts(oi(1), oi(2), oi(3)) = nb_sub_pts(oi(1), oi(2), oi(3))+1
                sub_pts(oi(1), oi(2), oi(3), nb_sub_pts(oi(1), oi(2), oi(3))) = pt_id
            end do
!
            do i = 1, nb_div
                do j = 1, nb_div
                    do k = 1, kmax
                        call this%children(i, j, k)%split(nb_pts, coordinates, dim, max_nb_pts, &
                                                          nb_sub_pts(i, j, k), sub_pts(i, j, k, :))
                    end do
                end do
            end do
!
            deallocate (sub_pts)
        end if
!
    end subroutine split_node
!
! ==================================================================================================
!
    subroutine get_pts_around(this, coor, distance, max_nb_pts, nb_pts, list_pts)
        class(Octree), intent(inout) :: this
        integer(kind=8), intent(in) :: max_nb_pts
        integer(kind=8), intent(inout) :: nb_pts, list_pts(max_nb_pts)
        real(kind=8), intent(in) :: coor(3), distance
!
        integer(kind=8) :: i_pt, oi(3), ocx, ocy, ocz, ocz_max
        real(kind=8) :: diff(3)
        type(OctreeNode) :: on
        aster_logical :: find
!
        nb_pts = 0
        on = this%node
        ocz_max = nb_div
        if (this%dim < 3) then
            ocz_max = 1
        end if
!
        find = ASTER_FALSE
        do while (.not. find)
            call on%find_children(coor, oi)
!
            if (on%children(oi(1), oi(2), oi(3))%nb_pts <= this%max_pts_by_node) then
                find = ASTER_TRUE
                if (on%children(oi(1), oi(2), oi(3))%is_strictly_inside(coor, distance)) then
                    do i_pt = 1, on%children(oi(1), oi(2), oi(3))%nb_pts
                        diff = coor-on%children(oi(1), oi(2), oi(3))%coordinates(:, i_pt)
                        if (norm2(diff) <= distance) then
                            nb_pts = nb_pts+1
                            list_pts(nb_pts) = on%children(oi(1), oi(2), oi(3))%list_pts(i_pt)
                        end if
                    end do
                else
                    do ocx = max(1, oi(1)-1), min(nb_div, oi(1)+1)
                        do ocy = max(1, oi(2)-1), min(nb_div, oi(2)+1)
                            do ocz = max(1, oi(3)-1), min(ocz_max, oi(3)+1)
                                do i_pt = 1, on%children(ocx, ocy, ocz)%nb_pts
                                    diff = coor-on%children(ocx, ocy, ocz)%coordinates(:, i_pt)
                                    if (norm2(diff) <= distance) then
                                        nb_pts = nb_pts+1
                                        list_pts(nb_pts) = on%children(ocx, ocy, ocz)%list_pts(i_pt)
                                    end if
                                end do
                            end do
                        end do
                    end do
                end if
            else
                on = on%children(oi(1), oi(2), oi(3))
            end if
        end do
!
        ASSERT(nb_pts <= max_nb_pts)
!
    end subroutine get_pts_around
!
end module
