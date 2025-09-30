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
subroutine export_octree(mesh_in, mesh_out, nb_max_pt, nb_max_level)
!
    use crea_maillage_module
    use octree_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/addGroupElem.h"
#include "asterfort/addGrpMa.h"
#include "asterfort/assert.h"
#include "asterfort/cargeo.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "MeshTypes_type.h"
!
    character(len=8), intent(in) :: mesh_in, mesh_out
    integer(kind=8), intent(in) :: nb_max_pt, nb_max_level
!
! ----------------------------------------------------------------------
!         TRANSFORMATION DES MAILLES POUR HHO
! ----------------------------------------------------------------------
! IN        mesh_in   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  mesh_out  K8  NOM DU MAILLAGE TRANSFORME
! IN        NBMA    I  NOMBRE DE MAILLES A TRAITER
! IN        LIMA    I  NUMERO ET TYPE DES MAILLES A TRAITER
! ----------------------------------------------------------------------
!
    type(Mmesh) :: mesh, mesh_oct
    type(Octree) :: oct
    integer(kind=8) :: i_node, node_id, i_level, i_cell, nb_cells, nb_cells_level
    integer(kind=8), allocatable :: cells_level(:), list_cells(:)
    real(kind=8), allocatable :: coordinates(:, :)
    character(len=24) :: name
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Create new mesh
!
    call mesh%init(mesh_in, convert_max=ASTER_FALSE)
!
    allocate (coordinates(3, mesh%nb_nodes))
    node_id = 0
    do i_node = 1, mesh%nb_total_nodes
        if (mesh%nodes(i_node)%keep) then
            node_id = node_id+1
            coordinates(:, node_id) = mesh%nodes(i_node)%coor
        end if
    end do
    ASSERT(node_id == mesh%nb_nodes)
!
! - Compute octree
!
    call oct%init(to_aster_int(mesh%nb_nodes), coordinates, nb_max_pt, nb_max_level)
    deallocate (coordinates)
    call mesh%clean()
!
! - create octree mesh
!
    call mesh_oct%converter%init()
    mesh_oct%dim_mesh = 3_4
    mesh_oct%max_nodes = mesh%nb_total_nodes
    mesh_oct%max_cells = 2_4*mesh%nb_total_nodes
    mesh_oct%isHPC = ASTER_FALSE
    mesh_oct%convert_max = ASTER_FALSE
!
    allocate (mesh_oct%nodes(mesh_oct%max_nodes))
    allocate (mesh_oct%renumber_nodes(mesh_oct%max_nodes))
    allocate (mesh_oct%cells(mesh_oct%max_cells))
    allocate (mesh_oct%renumber_cells(mesh_oct%max_cells))
!
! - Create mesh
!
    nb_cells = 0
    allocate (cells_level(8**min(8, nb_max_level)))
    call export_node(mesh_oct, oct%node, nb_cells, cells_level)
!
! - Fix cells
!
    call mesh_oct%fix(ASTER_FALSE, ASTER_FALSE, ASTER_FALSE, ASTER_TRUE, ASTER_FALSE, 1.d-10)
!
! - Copy mesh
!
    call mesh_oct%copy_mesh(mesh_out)
!
! - Add group
!
    call addGroupElem(mesh_out, nb_max_level+1)
    do i_level = 0, nb_max_level
        nb_cells_level = 0
        do i_cell = 1, nb_cells
            if (cells_level(i_cell) == i_level) then
                nb_cells_level = nb_cells_level+1
            end if
        end do
!
        if (nb_cells_level > 0) then
            allocate (list_cells(nb_cells_level))
            nb_cells_level = 0
            do i_cell = 1, nb_cells
                if (cells_level(i_cell) == i_level) then
                    nb_cells_level = nb_cells_level+1
                    list_cells(nb_cells_level) = i_cell
                end if
            end do
!
            name = "level_"
            call codent(i_level, 'G', name(7:), "F")
            call addGrpMa(mesh_out, name, list_cells, nb_cells_level)
            deallocate (list_cells)
        end if
    end do
    deallocate (cells_level)
!
! - Cleaning
!
    call mesh_oct%clean()
!
! - Update parameters for modified mesh (bounding box and dimensions)
    call cargeo(mesh_out)
!
    call jedema()
!

contains
    recursive subroutine export_node(moct, ocn, nb_cells_, cells_level_)
        type(Mmesh), intent(inout) :: moct
        type(OctreeNode), intent(in) :: ocn
        integer(kind=8), intent(inout) :: nb_cells_
        integer(kind=8), allocatable, intent(inout) :: cells_level_(:)
!
        integer(kind=8) :: ocx, ocy, ocz, ocz_max, nodes(8)
        real(kind=8) :: coor(3)
!
        if (ocn%is_leaf()) then
!
! - Add nodes
!
            coor = [ocn%xmin(1), ocn%xmin(2), ocn%xmin(3)]
            nodes(1) = to_aster_int(moct%add_node(coor, 1_4))
            coor = [ocn%xmin(1)+ocn%dx(1), ocn%xmin(2), ocn%xmin(3)]
            nodes(2) = to_aster_int(moct%add_node(coor, 1_4))
            coor = [ocn%xmin(1)+ocn%dx(1), ocn%xmin(2)+ocn%dx(2), ocn%xmin(3)]
            nodes(3) = to_aster_int(moct%add_node(coor, 1_4))
            coor = [ocn%xmin(1), ocn%xmin(2)+ocn%dx(2), ocn%xmin(3)]
            nodes(4) = to_aster_int(moct%add_node(coor, 1_4))
            if (ocn%dim == 3) then
                coor = [ocn%xmin(1), ocn%xmin(2), ocn%xmin(3)+ocn%dx(3)]
                nodes(5) = to_aster_int(moct%add_node(coor, 1_4))
                coor = [ocn%xmin(1)+ocn%dx(1), ocn%xmin(2), ocn%xmin(3)+ocn%dx(3)]
                nodes(6) = to_aster_int(moct%add_node(coor, 1_4))
                coor = [ocn%xmin(1)+ocn%dx(1), ocn%xmin(2)+ocn%dx(2), ocn%xmin(3)+ocn%dx(3)]
                nodes(7) = to_aster_int(moct%add_node(coor, 1_4))
                coor = [ocn%xmin(1), ocn%xmin(2)+ocn%dx(2), ocn%xmin(3)+ocn%dx(3)]
                nodes(8) = to_aster_int(moct%add_node(coor, 1_4))
            end if
!
! - Add cells - not save sub_entity to save memory
!
            nb_cells_ = nb_cells+1
            ASSERT(nb_cells <= size(cells_level_))
            cells_level_(nb_cells) = ocn%level
            if (ocn%dim == 3) then
                call moct%add_cell(MT_HEXA8, nodes(1:8), ASTER_FALSE)
            else
                call moct%add_cell(MT_QUAD4, nodes(1:4), ASTER_FALSE)
            end if
        else
            ocz_max = size(ocn%children, 3)
            if (ocn%dim < 3) then
                ocz_max = 1
            end if
!
            do ocx = 1, size(ocn%children, 1)
                do ocy = 1, size(ocn%children, 2)
                    do ocz = 1, ocz_max
                        call export_node(moct, ocn%children(ocx, ocy, ocz), &
                                         nb_cells_, cells_level_)
                    end do
                end do
            end do
        end if
!
    end subroutine

end subroutine
