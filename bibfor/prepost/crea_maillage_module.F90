! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: nicolas.pignet at edf.fr
!
module crea_maillage_module
!
!
implicit none
!
private
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/codlet.h"
#include "asterfort/cpclma.h"
#include "asterfort/elrfno.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jeccta.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/sdmail.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! Module to create a mesh for HHO
!
! --------------------------------------------------------------------------------------------------
!
!
    type Mconverter
        aster_logical :: to_convert(MT_NTYMAX) = ASTER_FALSE
        integer :: convert_to(MT_NTYMAX) = 0
        integer :: convert_max(MT_NTYMAX) = 0
        character(len=8) :: name(MT_NTYMAX) = ' '
        character(len=8) :: short_name(MT_NTYMAX) = ' '
        integer :: cata_type(MT_NTYMAX) = 0
        integer :: map_type(MT_NTYMAX) = 0
        integer :: dim(MT_NTYMAX) = 0
        integer :: nno(MT_NTYMAX) = 0
        integer :: nnos(MT_NTYMAX) = 0
! ----- member functions
        contains
        procedure, public, pass :: init => init_conv
        procedure, public, pass :: add_conversion
    end type
!
    type Mcell
        integer :: type = 0
        integer :: dim = 0
        integer :: id = 0
        integer :: nodes(27) = 0
        character(len=8) :: name = ' '
    end type
!
    type Mvolume
        integer :: type = 0
        integer :: nodes(27) = 0
        integer :: nb_faces = 0, faces(6) = 0
        integer :: nb_edges = 0, edges(12) = 0
    end type
!
    type Mface
        integer :: type = 0
        integer :: nnos = 0
        integer :: nodes(9) = 0
        integer :: nnos_sort(4) = 0
        integer :: nb_edges = 0, edges(4) = 0
        integer :: nb_volumes = 0, volumes(2) = 0
    end type
!
    type Medge
        integer :: type = 0
        integer :: nodes(4) = 0
        integer :: nnos_sort(2) = 0
!        integer :: nb_faces = 0, faces(10) = 0, max_faces = 10
    end type
!
    type Mnode
        integer :: id
        real(kind=8) :: coor(3) = 0.d0
        aster_logical :: keep = ASTER_FALSE
        aster_logical :: orphelan = ASTER_FALSE
        character(len=8) :: name = ' '
! used to improve search of edges and faces
! it could be improved a lot to decrease memory consumption
        integer :: max_faces = 0, nb_faces = 0
        integer :: max_edges = 0, nb_edges = 0
        integer, allocatable :: faces(:)
        integer, allocatable :: edges(:)
    end type
!
    type Mmesh
        integer :: nb_nodes = 0, nb_edges = 0, nb_faces = 0, nb_volumes = 0, nb_cells = 0
        integer :: nb_total_nodes = 0
        integer :: max_nodes = 0, max_edges = 0, max_faces = 0, max_volumes = 0, max_cells = 0
        integer :: dim_mesh = 0
!
        type(Mnode), allocatable :: nodes(:)
        type(Medge), allocatable :: edges(:)
        type(Mface), allocatable :: faces(:)
        type(Mvolume), allocatable :: volumes(:)
        type(Mcell), allocatable :: cells(:)
        type(Mconverter) :: converter
!
        character(len=8) :: mesh_in = ' '
        character(len=19) :: connex_in = ' '
!
        character(len=8)  :: node_prefix = 'N'
        integer :: node_index = 1
!
        integer, pointer :: v_typema(:) => null()
        aster_logical :: debug = ASTER_FALSE
        integer :: info = 0
! ----- member functions
        contains
        procedure, public, pass :: init => init_mesh
        procedure, public, pass :: clean => clean_mesh
        procedure, public, pass :: copy_mesh
        procedure, public, pass :: convert_cells
        procedure, public, pass :: add_cell
        procedure, public, pass :: check_mesh
        procedure, public, pass :: create_joints
        procedure, private, pass :: add_volume
        procedure, private, pass :: add_face
        procedure, private, pass :: add_edge
        procedure, private, pass :: add_node
        procedure, private, pass :: add_point1
        procedure, private, pass :: find_edge
        procedure, private, pass :: find_face
        procedure, private, pass :: convert_volume
        procedure, private, pass :: convert_face
        procedure, private, pass :: convert_edge
        procedure, private, pass :: barycenter
        procedure, private, pass :: update_nodes
        procedure, private, pass :: copy_group_no
        procedure, private, pass :: increase_memory
        procedure, private, pass :: numbering_nodes
    end type
!
!===================================================================================================
!
!===================================================================================================
!
    public :: Medge, Mface, Mcell, Mmesh, Mconverter
    private :: circ_perm, numbering_edge, numbering_face
contains
!
!===================================================================================================
!
!===================================================================================================
    subroutine add_conversion(this, from_type, to_type)
!
        implicit none
!
        class(Mconverter), intent(inout) :: this
        character(len=8), intent(in) :: from_type, to_type
!
        integer :: from_i, to_i
!
        call jemarq()
!
        call jenonu(jexnom('&CATA.TM.NOMTM', from_type), from_i)
        call jenonu(jexnom('&CATA.TM.NOMTM', to_type), to_i)
        ASSERT(from_i > 0)
        ASSERT(to_i > 0)
!
        this%to_convert(this%map_type(from_i)) = ASTER_TRUE
        this%convert_to(this%map_type(from_i)) = this%map_type(to_i)
!
        call jedema()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
    subroutine find_edge(this, nnos_sort, find, edge_id)
!
        implicit none
!
            class(Mmesh), intent(in) :: this
            integer, intent(in) :: nnos_sort(2)
            aster_logical, intent(out) :: find
            integer, intent(out) :: edge_id
!
! Find edge_id of edges
!
            integer :: i_edge, nb_edges, edge_i
            aster_logical :: ok
!
            find = ASTER_FALSE
            edge_id = 0
!
            nb_edges = this%nodes(nnos_sort(1))%nb_edges
            if(nb_edges > 0) then
                do i_edge = 1, nb_edges
                    edge_i = this%nodes(nnos_sort(1))%edges(i_edge)
                    ok = ASTER_TRUE
                    if( nnos_sort(1) .ne. this%edges(edge_i)%nnos_sort(1)) then
                        ok = ASTER_FALSE
                    else
                        if( nnos_sort(2) .ne. this%edges(edge_i)%nnos_sort(2)) then
                            ok = ASTER_FALSE
                        end if
                    end if
                    if(ok) then
                        find = ASTER_TRUE
                        edge_id = edge_i
                        exit
                    end if
                end do
            end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
    subroutine find_face(this, nnos, nnos_sort, find, face_id)
!
        implicit none
!
            class(Mmesh), intent(in) ::this
            integer, intent(in) :: nnos
            integer, intent(in) :: nnos_sort(4)
            aster_logical, intent(out) :: find
            integer, intent(out) :: face_id
!
! Find face_id of face
!
            integer :: i_face, i_node, nb_faces, face_i
            aster_logical :: ok
!
            find = ASTER_FALSE
            face_id = 0
!
            nb_faces = this%nodes(nnos_sort(1))%nb_faces
            if(nb_faces > 0) then
                do i_face = 1, nb_faces
                    face_i = this%nodes(nnos_sort(1))%faces(i_face)
                    if(nnos == this%faces(face_i)%nnos) then
                        ok = ASTER_TRUE
                        do i_node = 1, nnos
                            if( nnos_sort(i_node) .ne. this%faces(face_i)%nnos_sort(i_node)) then
                                ok = ASTER_FALSE
                                exit
                            end if
                        end do
                    else
                        ok = ASTER_FALSE
                    end if
                    if(ok) then
                        find = ASTER_TRUE
                        face_id = face_i
                        exit
                    end if
                end do
            end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine circ_perm(nb_nodes, nodes)
!
! Performs circular permutation such that the first element is the smallest and
! the second element is the second smallest
!
        implicit none
!
            integer, intent(in) :: nb_nodes
            integer, intent(inout) :: nodes(1:nb_nodes)
!
            integer :: tmp(27), ind_min, i_node, val_min
!
            ASSERT(nb_nodes <= 27)

            tmp(1:nb_nodes) = nodes(1:nb_nodes)
!
            ind_min = 1
            val_min = nodes(1)
            do i_node = 2, nb_nodes
                if(nodes(i_node) == val_min) then
                    ASSERT(ASTER_FALSE)
                else if(nodes(i_node) < val_min) then
                    ind_min = i_node
                    val_min = nodes(i_node)
                end if
            end do
            ASSERT(ind_min >= 1 .and. ind_min <= nb_nodes)
!
            if(ind_min > 1) then
                do i_node = ind_min, nb_nodes
                    nodes(i_node - ind_min + 1) = tmp(i_node)
                end do
                do i_node = 1, ind_min-1
                    nodes(i_node + nb_nodes - ind_min + 1) = tmp(i_node)
                end do
            end if
!
            if(nodes(nb_nodes) > nodes(2)) then
                tmp(1:nb_nodes) = nodes(1:nb_nodes)

                do i_node = 2, nb_nodes
                    nodes(i_node) = tmp(nb_nodes - i_node + 2)
                end do
            end if
!
            ASSERT(nodes(1) == val_min)
!
            ! print*,"AVANT: ", tmp(1:nb_elem)
            ! print*,"APRES: ", elems(1:nb_elem)
    end subroutine
!
! ==================================================================================================
!
    subroutine init_conv(this)
!
        implicit none
!
            class(Mconverter), intent(inout) :: this
!
            integer :: i_type, type_nume
!
            call jemarq()
!
            this%name(MT_POI1)     = "POI1"
            this%name(MT_SEG2)     = "SEG2"
            this%name(MT_SEG3)     = "SEG3"
            this%name(MT_SEG4)     = "SEG4"
            this%name(MT_TRIA3)    = "TRIA3"
            this%name(MT_TRIA6)    = "TRIA6"
            this%name(MT_TRIA7)    = "TRIA7"
            this%name(MT_QUAD4)    = "QUAD4"
            this%name(MT_QUAD8)    = "QUAD8"
            this%name(MT_QUAD9)    = "QUAD9"
            this%name(MT_TETRA4)   = "TETRA4"
            this%name(MT_TETRA10)  = "TETRA10"
            this%name(MT_TETRA15)  = "TETRA15"
            this%name(MT_HEXA8)    = "HEXA8"
            this%name(MT_HEXA9)    = "HEXA9"
            this%name(MT_HEXA20)   = "HEXA20"
            this%name(MT_HEXA27)   = "HEXA27"
            this%name(MT_PYRAM5)   = "PYRAM5"
            this%name(MT_PYRAM13)  = "PYRAM13"
            this%name(MT_PYRAM19)  = "PYRAM19"
            this%name(MT_PENTA6)   = "PENTA6"
            this%name(MT_PENTA7)   = "PENTA7"
            this%name(MT_PENTA15)  = "PENTA15"
            this%name(MT_PENTA18)  = "PENTA18"
            this%name(MT_PENTA21)  = "PENTA21"
!
            this%short_name(MT_POI1)     = "PO1"
            this%short_name(MT_SEG2)     = "SE2"
            this%short_name(MT_SEG3)     = "SE3"
            this%short_name(MT_SEG4)     = "SE4"
            this%short_name(MT_TRIA3)    = "TR3"
            this%short_name(MT_TRIA6)    = "TR6"
            this%short_name(MT_TRIA7)    = "TR7"
            this%short_name(MT_QUAD4)    = "QU4"
            this%short_name(MT_QUAD8)    = "QU8"
            this%short_name(MT_QUAD9)    = "QU9"
            this%short_name(MT_TETRA4)   = "TE4"
            this%short_name(MT_TETRA10)  = "T10"
            this%short_name(MT_TETRA15)  = "T15"
            this%short_name(MT_HEXA8)    = "HE8"
            this%short_name(MT_HEXA9)    = "HE9"
            this%short_name(MT_HEXA20)   = "H20"
            this%short_name(MT_HEXA27)   = "H27"
            this%short_name(MT_PYRAM5)   = "PY5"
            this%short_name(MT_PYRAM13)  = "P13"
            this%short_name(MT_PYRAM19)  = "P19"
            this%short_name(MT_PENTA6)   = "PE6"
            this%short_name(MT_PENTA7)   = "PE7"
            this%short_name(MT_PENTA15)  = "P15"
            this%short_name(MT_PENTA18)  = "P18"
            this%short_name(MT_PENTA21)  = "P21"
!
            this%convert_max(MT_POI1)     = MT_POI1
            this%convert_max(MT_SEG2)     = MT_SEG3
            this%convert_max(MT_SEG3)     = MT_SEG3
            this%convert_max(MT_SEG4)     = MT_SEG4
            this%convert_max(MT_TRIA3)    = MT_TRIA7
            this%convert_max(MT_TRIA6)    = MT_TRIA7
            this%convert_max(MT_TRIA7)    = MT_TRIA7
            this%convert_max(MT_QUAD4)    = MT_QUAD9
            this%convert_max(MT_QUAD8)    = MT_QUAD9
            this%convert_max(MT_QUAD9)    = MT_QUAD9
            this%convert_max(MT_TETRA4)   = MT_TETRA15
            this%convert_max(MT_TETRA10)  = MT_TETRA15
            this%convert_max(MT_TETRA15)  = MT_TETRA15
            this%convert_max(MT_HEXA8)    = MT_HEXA27
            this%convert_max(MT_HEXA9)    = MT_HEXA27
            this%convert_max(MT_HEXA20)   = MT_HEXA27
            this%convert_max(MT_HEXA27)   = MT_HEXA27
            this%convert_max(MT_PYRAM5)   = MT_PYRAM19
            this%convert_max(MT_PYRAM13)  = MT_PYRAM19
            this%convert_max(MT_PYRAM19)  = MT_PYRAM19
            this%convert_max(MT_PENTA6)   = MT_PENTA21
            this%convert_max(MT_PENTA7)   = MT_PENTA21
            this%convert_max(MT_PENTA15)  = MT_PENTA21
            this%convert_max(MT_PENTA18)  = MT_PENTA21
            this%convert_max(MT_PENTA21)  = MT_PENTA21
!
            do i_type = 1, MT_NTYMAX
                if(this%name(i_type) .ne. ' ') then
                    call jenonu(jexnom('&CATA.TM.NOMTM', this%name(i_type)), type_nume)
                    this%cata_type(i_type) = type_nume
                    this%map_type(type_nume) = i_type
                    call elrfno(this%short_name(i_type), this%nno(i_type), this%nnos(i_type), &
                                this%dim(i_type))
                    this%convert_to(i_type) = i_type
                end if
            end do
!
            call jedema()
!
    end subroutine
!
! ==================================================================================================
!
    subroutine init_mesh(this, mesh_in, info)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            character(len=8), intent(in) :: mesh_in
            integer, optional :: info
! --------------------------------------------------------------------------------------------------
! The idea is to read a given mesh (mesh_in). Internally, all the cells are stored with all nodes
! possibles for a cell. The internal cells stored are POI1, SEG3, SEG4, TRIA7, QUAD9,
! TETRA15, HEXA27, PENTA21 and PYRAM19
! For an other cell type is only necessary to know which are the nodes to use
! --------------------------------------------------------------------------------------------------
            integer, pointer :: v_mesh_dime(:) => null()
            real(kind=8), pointer :: v_coor(:) => null()
            character(len=8) :: name
            integer :: nb_elem_mesh, nb_node_mesh, i_node
            integer :: i_cell, nno, node_id
            real(kind=8):: start, end
!
            call jemarq()
            call this%converter%init()
!
            if(present(info)) then
                this%info = info
            end if
!
            if(this%info >= 2) then
                print*, "Creating mesh..."
                call cpu_time(start)
            end if
!
            call jeveuo(mesh_in//'.DIME', 'L', vi = v_mesh_dime)
            this%dim_mesh = v_mesh_dime(6)
            nb_elem_mesh = v_mesh_dime(3)
            nb_node_mesh = v_mesh_dime(1)
            this%mesh_in = mesh_in
            this%connex_in = mesh_in//'.CONNEX'
!
            this%max_cells = nb_elem_mesh
            this%max_volumes = nb_elem_mesh
            this%max_faces = 6*nb_elem_mesh
            this%max_edges = 12*nb_elem_mesh
            this%max_nodes = 4*nb_node_mesh
!
            allocate(this%cells(this%max_cells))
            allocate(this%volumes(this%max_volumes))
            allocate(this%faces(this%max_faces))
            allocate(this%edges(this%max_edges))
            allocate(this%nodes(this%max_nodes))
!
            call jeveuo(mesh_in//'.TYPMAIL', 'L', vi=this%v_typema)
            call jeveuo(mesh_in// '.COORDO    .VALE', 'L', vr = v_coor)
!
! --- Fill mesh
!
            this%nb_total_nodes = this%nb_nodes
            do i_node = 1, nb_node_mesh
                call jenuno(jexnum(mesh_in//'.NOMNOE', i_node), name)
                node_id = this%add_node(v_coor(3*(i_node-1)+1:3*(i_node-1)+3), name)
                ASSERT(i_node == node_id)
                this%nodes(node_id)%orphelan = ASTER_TRUE
                this%nodes(node_id)%max_faces = 30
                allocate(this%nodes(node_id)%faces(this%nodes(node_id)%max_faces))
                this%nodes(node_id)%max_edges = 30
                allocate(this%nodes(node_id)%edges(this%nodes(node_id)%max_edges))
            end do
!
            do i_cell = 1, nb_elem_mesh
                call this%add_cell(i_cell)
            end do
!
! --- Search orphelan nodes - to keep at the end
            do i_cell = 1, this%nb_cells
                nno = this%converter%nno(this%cells(i_cell)%type)
                do i_node = 1, nno
                    node_id = this%cells(i_cell)%nodes(i_node)
                    this%nodes(node_id)%orphelan = ASTER_FALSE
                end do
            end do
!
            if(this%info >= 2) then
                call cpu_time(end)
                print*, "... in ", end-start, " seconds."
            end if
!
            call jedema()
!
    end subroutine
!
! ==================================================================================================
!
    subroutine clean_mesh(this)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
!
            integer :: i_node, nb_nodes
!
            if(this%info >= 2) then
                print*, "Cleaning objects..."
            end if
!
            nb_nodes = size(this%nodes)
            do i_node = 1, nb_nodes
                if( allocated(this%nodes(i_node)%faces) ) then
                    deallocate(this%nodes(i_node)%faces)
                end if
                if( allocated(this%nodes(i_node)%edges) ) then
                    deallocate(this%nodes(i_node)%edges)
                end if
            end do
!
            deallocate(this%nodes)
            deallocate(this%edges)
            deallocate(this%faces)
            deallocate(this%volumes)
            deallocate(this%cells)
!
    end subroutine
!
! ==================================================================================================
!
    subroutine numbering_edge(cell_type, nb_edge, edge_type, edge_loc)
!
        implicit none
!
        integer, intent(in) :: cell_type
        integer, intent(out) :: nb_edge, edge_type(12), edge_loc(3,12)
! --------------------------------------------------------------------------------------------------
! Get edge connectivity of a cell
!
! --------------------------------------------------------------------------------------------------
        nb_edge = 0
        edge_type = 0
        edge_loc = 0
!
        if(cell_type == MT_QUAD4 .or. cell_type == MT_QUAD8 .or. cell_type == MT_QUAD9) then
            nb_edge = 4
            edge_loc(1:2,1) = [1,2]
            edge_loc(1:2,2) = [2,3]
            edge_loc(1:2,3) = [3,4]
            edge_loc(1:2,4) = [4,1]
            if(cell_type == MT_QUAD4) then
                edge_type(1:nb_edge) = MT_SEG2
            else
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3,1:nb_edge) = [5,6,7,8]
            end if
        elseif(cell_type == MT_TRIA3 .or. cell_type == MT_TRIA6 .or. cell_type == MT_TRIA7) then
            nb_edge = 3
            edge_loc(1:2,1) = [1,2]
            edge_loc(1:2,2) = [2,3]
            edge_loc(1:2,3) = [3,1]
            if(cell_type == MT_TRIA3) then
                edge_type(1:nb_edge) = MT_SEG2
            else
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3,1:nb_edge) = [4,5,6]
            end if
        elseif(cell_type == MT_HEXA8 .or. cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
            nb_edge = 12
            edge_loc(1:2,1) = [1,2]
            edge_loc(1:2,2) = [2,3]
            edge_loc(1:2,3) = [3,4]
            edge_loc(1:2,4) = [1,4]
            edge_loc(1:2,5) = [1,5]
            edge_loc(1:2,6) = [2,6]
            edge_loc(1:2,7) = [3,7]
            edge_loc(1:2,8) = [4,8]
            edge_loc(1:2,9) = [5,6]
            edge_loc(1:2,10) = [6,7]
            edge_loc(1:2,11) = [7,8]
            edge_loc(1:2,12) = [5,8]
            if(cell_type == MT_HEXA8) then
                edge_type(1:nb_edge) = MT_SEG2
           elseif(cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3,1:nb_edge) = [9,10,11,12,13,14,15,16,17,18,19,20]
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif(cell_type == MT_TETRA4 .or. cell_type == MT_TETRA10 .or. &
                cell_type == MT_TETRA15) then
            nb_edge = 6
            edge_loc(1:2,1) = [1,2]
            edge_loc(1:2,2) = [2,3]
            edge_loc(1:2,3) = [1,3]
            edge_loc(1:2,4) = [1,4]
            edge_loc(1:2,5) = [2,4]
            edge_loc(1:2,6) = [3,4]
            if(cell_type == MT_TETRA4) then
                edge_type(1:nb_edge) = MT_SEG2
            elseif(cell_type == MT_TETRA10 .or. cell_type == MT_TETRA15) then
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3,1:nb_edge) = [5,6,7,8,9,10]
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif(cell_type == MT_PENTA6 .or. cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 &
                .or. cell_type == MT_PENTA21) then
            nb_edge = 9
            edge_loc(1:2,1) = [1,2]
            edge_loc(1:2,2) = [2,3]
            edge_loc(1:2,3) = [1,3]
            edge_loc(1:2,4) = [1,4]
            edge_loc(1:2,5) = [2,5]
            edge_loc(1:2,6) = [3,6]
            edge_loc(1:2,7) = [4,5]
            edge_loc(1:2,8) = [5,6]
            edge_loc(1:2,9) = [4,6]
            if(cell_type == MT_PENTA6) then
                edge_type(1:nb_edge) = MT_SEG2
            elseif(cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 .or. &
                cell_type == MT_PENTA21) then
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3,1:nb_edge) = [7,8,9,10,11,12,13,14,15]
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif(cell_type == MT_PYRAM5 .or. cell_type == MT_PYRAM13 .or. &
                cell_type == MT_PYRAM19) then
            nb_edge = 8
            edge_loc(1:2,1) = [1,2]
            edge_loc(1:2,2) = [2,3]
            edge_loc(1:2,3) = [3,4]
            edge_loc(1:2,4) = [1,4]
            edge_loc(1:2,5) = [1,5]
            edge_loc(1:2,6) = [2,5]
            edge_loc(1:2,7) = [3,5]
            edge_loc(1:2,8) = [4,5]
            if(cell_type == MT_PYRAM5) then
                edge_type(1:nb_edge) = MT_SEG2
            elseif(cell_type == MT_PYRAM13 .or. cell_type == MT_PYRAM19) then
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3,1:nb_edge) = [6,7,8,9,10,11,12,13]
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

    end subroutine
!
! ==================================================================================================
!
    subroutine numbering_face(cell_type, nb_face, face_type, face_loc)
!
        implicit none
!
        integer, intent(in) :: cell_type
        integer, intent(out) :: nb_face, face_type(6), face_loc(9,6)
! ---------------------------------------------------------------------------------
! Get face connectivity of a cell
!
! ---------------------------------------------------------------------------------
        nb_face = 0
        face_type = 0
        face_loc = 0
!
        if(cell_type == MT_HEXA8 .or. cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
            nb_face = 6
            face_loc(1:4,1) = [1,4,3,2]
            face_loc(1:4,2) = [1,2,6,5]
            face_loc(1:4,3) = [2,3,7,6]
            face_loc(1:4,4) = [3,4,8,7]
            face_loc(1:4,5) = [1,5,8,4]
            face_loc(1:4,6) = [5,6,7,8]
            if(cell_type == MT_HEXA8) then
                face_type(1:6) = MT_QUAD4
           elseif(cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
                face_type(1:6) = MT_QUAD8
                face_loc(5:8,1) = [12,11,10,9]
                face_loc(5:8,2) = [9,14,17,13]
                face_loc(5:8,3) = [10,15,18,14]
                face_loc(5:8,4) = [11,16,19,15]
                face_loc(5:8,5) = [13,20,16,12]
                face_loc(5:8,6) = [17,18,19,20]
                if(cell_type == MT_HEXA27) then
                    face_type(1:6) = MT_QUAD9
                    face_loc(9,1:6) = [21,22,23,24,25,26]
                endif
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif(cell_type == MT_TETRA4 .or. cell_type == MT_TETRA10 &
                .or. cell_type == MT_TETRA15) then
            nb_face = 4
            face_loc(1:3,1) = [1,3,2]
            face_loc(1:3,2) = [1,2,4]
            face_loc(1:3,3) = [1,4,3]
            face_loc(1:3,4) = [2,3,4]
            if(cell_type == MT_TETRA4) then
                face_type(1:4) = MT_TRIA3
            elseif(cell_type == MT_TETRA10 .or. cell_type == MT_TETRA15) then
                face_type(1:4) = MT_TRIA6
                face_loc(4:6,1) = [7,6,5]
                face_loc(4:6,2) = [1,9,8]
                face_loc(4:6,3) = [8,10,7]
                face_loc(4:6,4) = [6,10,9]
                if(cell_type == MT_TETRA15) then
                    face_type(1:4) = MT_TRIA7
                    face_loc(7,1:4) = [11,12,13,14]
                endif
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif(cell_type == MT_PENTA6 .or. cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 &
            .or. cell_type == MT_PENTA21) then
            nb_face = 5
            face_loc(1:4,1) = [1,2,5,4]
            face_loc(1:4,2) = [2,3,6,5]
            face_loc(1:4,3) = [1,4,6,3]
            face_loc(1:3,4) = [1,3,2]
            face_loc(1:3,5) = [4,5,6]
            if(cell_type == MT_PENTA6) then
                face_type(1:3) = MT_QUAD4
                face_type(4:5) = MT_TRIA3
            elseif(cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 .or. &
                    cell_type == MT_PENTA21) then
                face_type(1:3) = MT_QUAD8
                face_type(4:5) = MT_TRIA6
                face_loc(5:8,1) = [7,11,13,10]
                face_loc(5:8,2) = [8,12,14,11]
                face_loc(5:8,3) = [10,15,12,9]
                face_loc(4:6,4) = [9,8,7]
                face_loc(4:6,5) = [13,14,15]
                if(cell_type == MT_PENTA18 .or. cell_type == MT_PENTA21) then
                    face_type(1:3) = MT_QUAD9
                    face_loc(9,1:3) = [16,17,18]
                    if(cell_type == MT_PENTA21) then
                        face_type(4:5) = MT_TRIA7
                        face_loc(7,4:5) = [19,20]
                    endif
                endif
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif(cell_type == MT_PYRAM5 .or. cell_type == MT_PYRAM13 .or. &
                cell_type == MT_PYRAM19) then
            nb_face = 5
            face_loc(1:4,1) = [1,2,3,4]
            face_loc(1:3,2) = [1,2,5]
            face_loc(1:3,3) = [2,3,5]
            face_loc(1:3,4) = [3,4,5]
            face_loc(1:3,5) = [4,1,5]
            if(cell_type == MT_PYRAM5) then
                face_type(1) = MT_QUAD4
                face_type(2:5) = MT_TRIA3
            elseif(cell_type == MT_PYRAM13 .or. cell_type == MT_PYRAM19) then
                face_type(1) = MT_QUAD8
                face_type(2:5) = MT_TRIA6
                face_loc(5:8,1) = [6,7,8,9]
                face_loc(4:6,2) = [6,11,10]
                face_loc(4:6,3) = [7,12,11]
                face_loc(4:6,4) = [8,13,12]
                face_loc(4:6,5) = [9,10,13]
                if(cell_type == MT_PYRAM19) then
                    face_type(1) = MT_QUAD9
                    face_type(2:5) = MT_TRIA7
                    face_loc(9,1) = 14
                    face_loc(7,2:5) = [15,16,17,18]
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine numbering_nodes(this, cell_type, nodes_loc)
!
        implicit none
!
        class(Mmesh), intent(in) :: this
        integer, intent(in) :: cell_type
        integer, intent(out) :: nodes_loc(27)
! -----------------------------------------------------------------------------------
! Add special treatment
        integer :: i_node, nno
!
        nodes_loc = 0
        nno = this%converter%nno(cell_type)
!
        do i_node = 1, nno
            nodes_loc(i_node) = i_node
        end do
        if (cell_type .eq. MT_HEXA9) then
            nodes_loc(1:8) = [1, 2, 3, 4, 5, 6, 7, 8]
            nodes_loc(9)   = 27
        endif
        if (cell_type .eq. MT_PENTA7) then
            nodes_loc(1:6) = [1, 2, 3, 4, 5, 6]
            nodes_loc(7)   = 21
        endif
!
    end subroutine
!
! ==================================================================================================
!
    subroutine add_cell(this, cell_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: cell_id
!
            integer :: cell_type, cell_dim, cell_nodes(27), nb_nodes, cell_index
            integer, pointer :: v_connex(:) => null()
!
            ASSERT(this%nb_cells < this%max_cells)
            this%nb_cells = this%nb_cells + 1
            cell_type = this%converter%map_type(this%v_typema(cell_id))
            cell_dim = this%converter%dim(cell_type)
!
            call jeveuo(jexnum(this%connex_in, cell_id), 'L',  vi = v_connex)
            nb_nodes = this%converter%nno(cell_type)
            cell_nodes(1:nb_nodes) = v_connex(1:nb_nodes)
!
            if(this%debug) then
                print*, "Cell ", this%nb_cells, ": ", cell_type, &
                this%converter%name(cell_type), cell_dim
            end if
!
            if(cell_dim == 3) then
                cell_index = this%add_volume(cell_type, cell_nodes)
            elseif(cell_dim == 2) then
                cell_index = this%add_face(cell_type, cell_nodes)
            elseif(cell_dim == 1) then
                cell_index = this%add_edge(cell_type, cell_nodes)
            elseif(cell_dim == 0) then
                cell_index = this%add_point1(cell_type, cell_nodes)
            else
                ASSERT(ASTER_FALSE)
            end if
!
            this%cells(this%nb_cells)%type = cell_type
            this%cells(this%nb_cells)%dim = cell_dim
            this%cells(this%nb_cells)%id = cell_index
            this%cells(this%nb_cells)%nodes(1:nb_nodes) = cell_nodes(1:nb_nodes)
            call jenuno(jexnum(this%mesh_in//'.NOMMAI', cell_id), this%cells(this%nb_cells)%name)
!
    end subroutine
!
! ==================================================================================================
!
    function add_volume(this, type, nodes) result(volume_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: type, nodes(27)
            integer :: volume_id
! ----------------------------------------------------------------------
            integer :: nno, i_face
            integer :: nb_edges, edge_type(12), edge_loc(3,12), edge_id, edge_nno
            integer :: nb_faces, face_type(6), face_loc(9,6), i_node, i_edge
            integer :: face_nno, face_nodes(27), face_id, edge_nodes(27)
            aster_logical :: find
!
            ASSERT(this%converter%dim(type) == 3)
            nno = this%converter%nno(type)
!
            this%nb_volumes = this%nb_volumes + 1
            if(this%nb_volumes > this%max_volumes) then
                call this%increase_memory("VOLUMES ", 2*this%max_volumes)
            end if
            ASSERT(this%nb_volumes <= this%max_volumes)
            volume_id = this%nb_volumes
            this%volumes(volume_id)%type = type
            this%volumes(volume_id)%nodes(1:nno) = nodes(1:nno)
! --- create edges
            call numbering_edge(type, nb_edges, edge_type, edge_loc)
            this%volumes(volume_id)%nb_edges = nb_edges
            do i_edge = 1, nb_edges
                edge_nno = this%converter%nno(edge_type(i_edge))
                do i_node = 1, edge_nno
                    edge_nodes(i_node) = nodes(edge_loc(i_node, i_edge))
                end do
                edge_id = this%add_edge(edge_type(i_edge), edge_nodes, volume_id)
                this%volumes(volume_id)%edges(i_edge) = edge_id
            end do
! --- create faces
            call numbering_face(type, nb_faces, face_type, face_loc)
            this%volumes(volume_id)%nb_faces = nb_faces
            do i_face = 1, nb_faces
                face_nno = this%converter%nno(face_type(i_face))
                do i_node = 1, face_nno
                    face_nodes(i_node) = nodes(face_loc(i_node, i_face))
                end do
                face_id = this%add_face(face_type(i_face), face_nodes, volume_id)
                this%volumes(volume_id)%faces(i_face) = face_id
            end do
!
            call this%convert_volume(volume_id)
!
            if(this%debug) then
                print*, "Add volume: ", volume_id
                print*, "- Find: ", find
                print*, "- Type: ", this%converter%name(this%volumes(volume_id)%type), &
                                    "(",this%volumes(volume_id)%type, ")"
                print*, "- Nodes: ", this%volumes(volume_id)%nodes
                print*, "- Edges: ", this%volumes(volume_id)%edges
                print*, "- Faces: ", this%volumes(volume_id)%faces
            end if
!
    end function
!
! ==================================================================================================
!
    function add_face(this, type, nodes, volume_id) result(face_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: type, nodes(27)
            integer, intent(in), optional :: volume_id
            integer :: face_id
! ----------------------------------------------------------------------
            integer :: nno, nnos, nnos_sort(4), i_edge, edge_id
            integer :: nb_edge, edge_type(12), edge_loc(3,12), i_node
            integer :: edge_nno, edge_nodes(27), old_size
            integer, allocatable :: new_faces(:)
            aster_logical :: find
!
            ASSERT(this%converter%dim(type) == 2)
            nno = this%converter%nno(type)
            nnos = this%converter%nnos(type)
            nnos_sort(1:nnos) = nodes(1:nnos)
!
            call circ_perm(nnos, nnos_sort)
!
            call this%find_face(nnos, nnos_sort, find, face_id)
!
            if(.not. find) then
                this%nb_faces = this%nb_faces + 1
                if(this%nb_faces > this%max_faces) then
                    call this%increase_memory("FACES   ", 2*this%max_faces)
                end if
                ASSERT(this%nb_faces <= this%max_faces)
                face_id = this%nb_faces
                this%faces(face_id)%type = type
                this%faces(face_id)%nodes(1:nno) = nodes(1:nno)
                this%faces(face_id)%nnos = nnos
                this%faces(face_id)%nnos_sort(1:nnos) = nnos_sort(1:nnos)
! --- create edges
                call numbering_edge(type, nb_edge, edge_type, edge_loc)
                this%faces(face_id)%nb_edges = nb_edge
                do i_edge = 1, nb_edge
                    edge_nno = this%converter%nno(edge_type(i_edge))
                    do i_node = 1, edge_nno
                        edge_nodes(i_node) = nodes(edge_loc(i_node, i_edge))
                    end do
                    edge_id = this%add_edge(edge_type(i_edge), edge_nodes, face_id)
                    this%faces(face_id)%edges(i_edge) = edge_id
                end do
!
                if(this%nodes(nnos_sort(1))%nb_faces >= this%nodes(nnos_sort(1))%max_faces) then
                    old_size = this%nodes(nnos_sort(1))%max_faces
                    allocate(new_faces(old_size))
                    new_faces(1:old_size) = this%nodes(nnos_sort(1))%faces(1:old_size)
                    deallocate(this%nodes(nnos_sort(1))%faces)
                    allocate(this%nodes(nnos_sort(1))%faces(2*old_size))
                    this%nodes(nnos_sort(1))%faces(1:old_size) = new_faces(1:old_size)
                    deallocate(new_faces)
                    this%nodes(nnos_sort(1))%max_faces = 2 * old_size
                end if
                this%nodes(nnos_sort(1))%nb_faces = this%nodes(nnos_sort(1))%nb_faces + 1
                this%nodes(nnos_sort(1))%faces(this%nodes(nnos_sort(1))%nb_faces) = face_id
!
                call this%convert_face(face_id)
            end if
!
            if(present(volume_id)) then
                this%faces(face_id)%nb_volumes = this%faces(face_id)%nb_volumes + 1
                this%faces(face_id)%volumes(this%faces(face_id)%nb_volumes) = volume_id
            end if
!
            if(this%debug) then
                print*, "Add face: ", face_id
                print*, "- Find: ", find
                print*, "- Type: ", this%converter%name(this%faces(face_id)%type), &
                                    "(",this%faces(face_id)%type, ")"
                print*, "- Nodes: ", this%faces(face_id)%nodes
                print*, "- NNOS: ", this%faces(face_id)%nnos_sort
                print*, "- Edges: ", this%faces(face_id)%edges
                print*, "- Volumes: ", this%faces(face_id)%volumes
            end if
!
    end function
!
! ==================================================================================================
!
    function add_edge(this, type, nodes, face_id) result(edge_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: type, nodes(27)
            integer, intent(in), optional :: face_id
            integer :: edge_id
! ----------------------------------------------------------------------
            integer :: nno, nnos_sort(2), old_size
            integer, allocatable :: new_edges(:)
            aster_logical :: find
!
            ASSERT(this%converter%dim(type) == 1)
            nno = this%converter%nno(type)
            nnos_sort(1:2) = nodes(1:2)
!
            call circ_perm(2, nnos_sort)
!
            call this%find_edge(nnos_sort, find, edge_id)
!
            if(.not. find) then
                this%nb_edges = this%nb_edges + 1
                if(this%nb_edges > this%max_edges) then
                    call this%increase_memory("EDGES   ", 2*this%max_edges)
                end if
                ASSERT(this%nb_edges <= this%max_edges)
                edge_id = this%nb_edges
                this%edges(edge_id)%type = type
                this%edges(edge_id)%nodes(1:nno) = nodes(1:nno)
                this%edges(edge_id)%nnos_sort = nnos_sort
!
                if(this%nodes(nnos_sort(1))%nb_edges >= this%nodes(nnos_sort(1))%max_edges) then
                    old_size = this%nodes(nnos_sort(1))%max_edges
                    allocate(new_edges(old_size))
                    new_edges(1:old_size) = this%nodes(nnos_sort(1))%edges(1:old_size)
                    deallocate(this%nodes(nnos_sort(1))%edges)
                    allocate(this%nodes(nnos_sort(1))%edges(2*old_size))
                    this%nodes(nnos_sort(1))%edges(1:old_size) = new_edges(1:old_size)
                    deallocate(new_edges)
                    this%nodes(nnos_sort(1))%max_edges = 2 * old_size
                end if
                this%nodes(nnos_sort(1))%nb_edges = this%nodes(nnos_sort(1))%nb_edges + 1
                this%nodes(nnos_sort(1))%edges(this%nodes(nnos_sort(1))%nb_edges) = edge_id
!
                call this%convert_edge(edge_id)
            end if
!
            if(present(face_id)) then
                ! this%edges(face_id)%nb_faces = this%edges(face_id)%nb_faces + 1
                ! ASSERT(this%edges(face_id)%nb_faces <= this%edges(face_id)%max_faces)
                ! this%edges(face_id)%faces(this%edges(face_id)%nb_faces) = face_id
            end if
!
            if(this%debug) then
                print*, "Add edge: ", edge_id
                print*, "- Find: ", find
                print*, "- Type: ", this%converter%name(this%edges(edge_id)%type), &
                                    "(",this%edges(edge_id)%type, ")"
                print*, "- Nodes: ", this%edges(edge_id)%nodes
                print*, "- NNOS: ", this%edges(edge_id)%nnos_sort
            end if
!
    end function
!
! ==================================================================================================
!
    function add_point1(this, type, nodes) result(point_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: type, nodes(27)
            integer :: point_id
! ----------------------------------------------------------------------
!
            ASSERT(this%converter%dim(type) == 0)
            ASSERT(this%converter%nno(type) == 1)
            point_id = 1
!
            if(this%debug) then
                print*, "Add point1: "
                print*, "- Type: ", this%converter%name(type), "(", type, ")"
                print*, "- Nodes: ", nodes(1)
                print*, "- NNOS: ", nodes(1)
            end if
!
    end function
!
! ==================================================================================================
!
    function add_node(this, coor, name) result(node_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            real(kind=8), intent(in) :: coor(3)
            character(len=8), intent(in), optional :: name
            integer :: node_id
! ----------------------------------------------------------------------
!
            this%nb_nodes = this%nb_nodes + 1
            this%nb_total_nodes = this%nb_total_nodes + 1
            if(this%nb_nodes > this%max_nodes) then
                call this%increase_memory("NODES   ", 2*this%max_nodes)
            end if
            ASSERT(this%nb_nodes <= this%max_nodes)
            node_id = this%nb_nodes
            this%nodes(node_id)%id = node_id
            this%nodes(node_id)%keep = ASTER_TRUE
            this%nodes(node_id)%coor(1:3) = coor
            if(present(name)) then
                this%nodes(node_id)%name = name
            else
                this%nodes(node_id)%name = 'XXXXXXXX'
            end if
!
    end function
!
! ==================================================================================================
!
    subroutine copy_mesh(this, mesh_out)
!
        implicit none
!
            class(Mmesh), intent(in) :: this
            character(len=8), intent(in) :: mesh_out
! ------------------------------------------------------------------
            character(len=24) :: nommai, nomnoe, cooval, coodsc, cooref, grpnoe
            character(len=24) :: gpptnn, grpmai, gpptnm, connex, titre, typmai, adapma
            character(len=4) :: dimesp
            integer :: i_node, nno, i_cell, ntgeo, nbnoma, node_id
            real(kind=8):: start, end
            real(kind=8), pointer :: v_coor(:) => null()
            integer, pointer :: v_int(:) => null()
            integer, pointer :: v_connex(:) => null()
            character(len=24), pointer :: v_k24(:) => null()
!
            call jemarq()
!
            call this%check_mesh()
!
            if(this%info >= 2) then
                print*, "Copying mesh..."
                call cpu_time(start)
            end if
!
            call sdmail(mesh_out, nommai, nomnoe, cooval, coodsc,&
                        cooref, grpnoe, gpptnn, grpmai, gpptnm,&
                        connex, titre, typmai, adapma)
!
! --- Create nodes
!
! ------ Set names
            call jecreo(nomnoe, 'G N K8')
            call jeecra(nomnoe, 'NOMMAX', this%nb_nodes)
            do i_node = 1 , this%nb_total_nodes
                if(this%nodes(i_node)%keep) then
                    call jecroc(jexnom(nomnoe, this%nodes(i_node)%name))
                end if
            end do
! ------ Copy coordinates
            call wkvect(cooval, 'G V R', this%nb_nodes*3, vr=v_coor)
            call codent(this%dim_mesh, 'G', dimesp)
            call jeecra(cooval, 'DOCU', cval=dimesp)
            node_id = 0
            do i_node = 1, this%nb_total_nodes
                if(this%nodes(i_node)%keep) then
                    node_id = node_id + 1
                    ASSERT(node_id == this%nodes(i_node)%id)
                    v_coor(3*(node_id-1)+1:3*(node_id-1)+3) = this%nodes(i_node)%coor(1:3)
                end if
            end do
! ------ Type of GEOM_R field
            call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), ntgeo)
            call wkvect(coodsc, 'G V I', 3, vi=v_int)
            call jeecra(coodsc, 'DOCU', 0, 'CHNO')
            v_int(1) = ntgeo
            v_int(2) = -3
            v_int(3) = 14
!
            call wkvect(cooref, 'G V K24', 4, vk24=v_k24)
            v_k24(1) = mesh_out
!
! --- Create cells
!
! ------ Set names
            call jecreo(nommai, 'G N K8')
            call jeecra(nommai, 'NOMMAX', this%nb_cells)
            do i_cell = 1, this%nb_cells
                call jecroc(jexnom(nommai, this%cells(i_cell)%name))
            end do
! ------ Count total number of nodes (repeated)
            nbnoma = 0
            do i_cell = 1, this%nb_cells
                nbnoma = nbnoma + this%converter%nno(this%cells(i_cell)%type)
            end do
! ------ Create connectivity
            call wkvect(typmai, 'G V I', this%nb_cells, vi=v_int)
            call jecrec(connex, 'G V I', 'NU', 'CONTIG', 'VARIABLE', this%nb_cells)
            call jeecra(connex, 'LONT', nbnoma)
!
            do i_cell = 1, this%nb_cells
                v_int(i_cell) = this%converter%cata_type(this%cells(i_cell)%type)
                nno = this%converter%nno(this%cells(i_cell)%type)
                call jeecra(jexnum(connex, i_cell), 'LONMAX', nno)
                call jeveuo(jexnum(connex, i_cell), 'E', vi=v_connex)
                do i_node = 1, nno
                    v_connex(i_node) = this%nodes(this%cells(i_cell)%nodes(i_node))%id
                end do
            end do
!
            call jeccta(connex)
!
! --- Create groups
!
            call this%copy_group_no(grpnoe, gpptnn)
            call cpclma(this%mesh_in, mesh_out, 'GROUPEMA', 'G')
!
! --- Create .DIME
!
            call wkvect(mesh_out//'.DIME', 'G V I', 6, vi=v_int)
            v_int(1)= this%nb_nodes
            v_int(3)= this%nb_cells
            v_int(6)= this%dim_mesh
!
            if(this%info >= 2) then
                call cpu_time(end)
                print*, "... in ", end-start, " seconds."
            end if
!
            call jedema()
!
    end subroutine
!
! ==================================================================================================
!
    subroutine convert_cells(this, nb_cells, list_cells, prefix, ndinit)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: nb_cells, list_cells(nb_cells)
            character(len=8), intent(in), optional :: prefix
            integer, intent(in), optional :: ndinit
! ------------------------------------------------------------------
            integer :: i_cell, cell_id, cell_dim, object_id, nno, cell_type
            integer :: i_node, nodes_loc(27)
            real(kind=8):: start, end
!
            if(this%info >= 2) then
                print*, "Converting cells..."
                call cpu_time(start)
            end if
!
            if(present(prefix)) then
                this%node_prefix = prefix
            end if
            if(present(ndinit)) then
                this%node_index = ndinit
            end if
!
            do i_cell = 1, nb_cells
                cell_id = list_cells(i_cell)
                cell_dim = this%cells(cell_id)%dim
                cell_type = this%cells(cell_id)%type
                object_id = this%cells(cell_id)%id
!
                if(this%debug) then
                    print*, "Convert ", cell_id, ": ", this%cells(cell_id)%type, &
                        this%converter%name(this%cells(cell_id)%type), cell_dim
                end if
!
                if(this%converter%to_convert(cell_type)) then
                    this%cells(cell_id)%type = this%converter%convert_to(cell_type)
                    nno = this%converter%nno(this%cells(cell_id)%type)
                    call this%numbering_nodes(this%cells(cell_id)%type, nodes_loc)
                    this%cells(cell_id)%nodes = 0
                    if(cell_dim == 3) then
                        do i_node = 1, nno
                            this%cells(cell_id)%nodes(i_node) = &
                                this%volumes(object_id)%nodes(nodes_loc(i_node))
                        end do
                    elseif(cell_dim == 2) then
                        do i_node = 1, nno
                            this%cells(cell_id)%nodes(i_node) = &
                                this%faces(object_id)%nodes(nodes_loc(i_node))
                        end do
                    elseif(cell_dim == 1) then
                        do i_node = 1, nno
                            this%cells(cell_id)%nodes(i_node) = &
                                this%edges(object_id)%nodes(nodes_loc(i_node))
                        end do
                    elseif(cell_dim == 0) then
                        ASSERT(ASTER_FALSE)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                end if
            end do
! --- Keep only necessary nodes
            call this%update_nodes()
!
            if(this%info >= 2) then
                call cpu_time(end)
                print*, "... in ", end-start, " seconds."
            end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine convert_volume(this, volu_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: volu_id
! ------------------------------------------------------------------
            integer :: volu_type, volu_type_end
            integer :: nno, nno_end, node_id, i_face, i_node, face_type, face_nno
            integer :: nb_edges, edge_type(12), edge_loc(3,12), face_id
            integer :: nb_faces, faces_type(6), face_loc(9,6), i_edge
!
            volu_type = this%volumes(volu_id)%type
            nno = this%converter%nno(volu_type)
            volu_type_end = this%converter%convert_max(volu_type)
            nno_end = this%converter%nno(volu_type_end)
!
            if(this%debug) then
                print*, "Volume: ", volu_id, volu_type, volu_type_end
            end if
!
            this%volumes(volu_id)%type = volu_type_end
!
            if(nno_end > nno) then
! --- Add nodes on edges
                if(volu_type == MT_HEXA8 .or. volu_type == MT_TETRA4 .or.  &
                    volu_type == MT_PYRAM5 .or. volu_type == MT_PENTA6) then
                    call numbering_edge(volu_type_end, nb_edges, edge_type, edge_loc)
                    do i_edge = 1, this%volumes(volu_id)%nb_edges
                        call this%convert_edge(this%volumes(volu_id)%edges(i_edge))
                        this%volumes(volu_id)%nodes(edge_loc(3, i_edge)) = &
                        this%edges(this%volumes(volu_id)%edges(i_edge))%nodes(3)
                    end do
                end if
! --- Add nodes on faces
                call numbering_face(volu_type_end, nb_faces, faces_type, face_loc)
                do i_face = 1, this%volumes(volu_id)%nb_faces
                    face_id = this%volumes(volu_id)%faces(i_face)
                    call this%convert_face(this%volumes(volu_id)%faces(i_face))
                    face_type = faces_type(i_face)
                    face_nno = this%converter%nno(face_type)
!
                    this%volumes(volu_id)%nodes(face_loc(face_nno, i_face)) = &
                                    this%faces(face_id)%nodes(face_nno)
                end do
! --- Add node at barycenter
                node_id = this%add_node(this%barycenter(nno_end-1, this%volumes(volu_id)%nodes))
                this%volumes(volu_id)%nodes(nno_end) = node_id
            else
                ASSERT(volu_type == volu_type_end)
            end if
!
            if(this%debug) then
                do i_node = 1, nno_end
                    ASSERT(this%volumes(volu_id)%nodes(i_node) > 0)
                end do
            end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine convert_face(this, face_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: face_id
! ------------------------------------------------------------------
            integer :: face_type, face_type_end
            integer :: nno, nno_end, node_id, i_edge, i_node
!
            face_type = this%faces(face_id)%type
            nno = this%converter%nno(face_type)
            face_type_end = this%converter%convert_max(face_type)
            nno_end = this%converter%nno(face_type_end)
!
            if(this%debug) then
                print*, "Face: ", face_id, face_type, face_type_end, this%faces(face_id)%nodes
            end if
!
            this%faces(face_id)%type = face_type_end
!
            if(nno_end > nno) then
! --- Add nodes at the middle of edges
                if(face_type == MT_TRIA3 .or. face_type == MT_QUAD4) then
                    do i_edge = 1, this%faces(face_id)%nb_edges
                        call this%convert_edge(this%faces(face_id)%edges(i_edge))
                        this%faces(face_id)%nodes(nno+i_edge) = &
                             this%edges(this%faces(face_id)%edges(i_edge))%nodes(3)
                    end do
                end if
! --- Add node at the barycenter
                node_id = this%add_node(this%barycenter(nno_end-1, this%faces(face_id)%nodes))
                this%faces(face_id)%nodes(nno_end) = node_id
            else
                ASSERT(face_type == face_type_end)
            end if
!
            if(this%debug) then
                do i_node = 1, nno_end
                    ASSERT(this%faces(face_id)%nodes(i_node) > 0)
                end do
            end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine convert_edge(this, edge_id)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            integer, intent(in) :: edge_id
! ------------------------------------------------------------------
            integer :: edge_type, edge_type_end
            integer :: nno, nno_end, node_id, i_node
!
            edge_type = this%edges(edge_id)%type
            nno = this%converter%nno(edge_type)
            edge_type_end = this%converter%convert_max(edge_type)
            nno_end = this%converter%nno(edge_type_end)
!
            if(this%debug) then
                print*, "Edge: ", edge_id, edge_type, edge_type_end, this%edges(edge_id)%nodes
            end if
!
            this%edges(edge_id)%type = edge_type_end
!
            if(nno_end > nno) then
                if(edge_type_end == MT_SEG3) then
                    node_id = this%add_node(this%barycenter(2, this%edges(edge_id)%nodes))
                    this%edges(edge_id)%nodes(3) = node_id
                else
                    ASSERT(ASTER_FALSE)
                    node_id = this%add_node([0.d0, 0.d0, 0.d0])
                    this%edges(edge_id)%nodes(3) = node_id
                    node_id = this%add_node([0.d0, 0.d0, 0.d0]  )
                    this%edges(edge_id)%nodes(4) = node_id
                end if
            else
                ASSERT(edge_type == edge_type_end)
            end if
!
            if(this%debug) then
                do i_node = 1, nno_end
                    ASSERT(this%edges(edge_id)%nodes(i_node) > 0)
                end do
            end if
!
    end subroutine
!
! ==================================================================================================
!
    function barycenter(this, nb_nodes, nodes) result(coor)
!
        implicit none
!
            class(Mmesh), intent(in) :: this
            integer, intent(in) :: nb_nodes, nodes(nb_nodes)
            real(kind=8) :: coor(3)
! ---------------------------------------------------------------------------------
            integer :: i_node
!
            coor = 0.d0
            do i_node = 1, nb_nodes
                coor(1:3) = coor(1:3) + this%nodes(nodes(i_node))%coor(1:3)
            end do
            coor(1:3) = coor(1:3) / real(nb_nodes, kind=8)
    end function
!
! ==================================================================================================
!
    subroutine update_nodes(this)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
! -----------------------------------------------------------------------
            integer :: i_node, i_cell, nno, node_id
            character(len=8) :: nume
!
            do i_node = 1, this%nb_total_nodes
                this%nodes(i_node)%keep = ASTER_FALSE
            end do
!
! --- Keep only nodes of cells
            do i_cell = 1, this%nb_cells
                nno = this%converter%nno(this%cells(i_cell)%type)
!
                if(this%debug) then
                    print*, "Cell: ", i_cell, this%cells(i_cell)%type, nno, &
                    this%cells(i_cell)%nodes(1:nno)
                end if
!
                do i_node = 1, nno
                    node_id = this%cells(i_cell)%nodes(i_node)
                    this%nodes(node_id)%keep = ASTER_TRUE
                end do
            end do
!
! --- Keep initial orphelan nodes
            do i_node = 1, this%nb_total_nodes
                if(this%nodes(i_node)%orphelan) then
                    this%nodes(i_node)%keep = ASTER_TRUE
                end if
            end do
!
! --- Renumbering and rename
            this%nb_nodes = 0
            do i_node = 1, this%nb_total_nodes
                if(this%nodes(i_node)%keep) then
                    this%nb_nodes = this%nb_nodes + 1
                    this%nodes(i_node)%id = this%nb_nodes
                    if(this%nodes(i_node)%name == "XXXXXXXX") then
                        if(this%node_index .ge. 10000000) then
                            call codlet(this%node_index, 'G', nume)
                        else
                            call codent(this%node_index, 'G', nume)
                        endif
                        this%nodes(i_node)%name = trim(this%node_prefix)//trim(nume)
                        this%node_index = this%node_index + 1
                    end if
                end if
            end do
!
            if(this%debug) then
                print*, "Update nodes: ", this%nb_nodes, this%nb_total_nodes
                do i_node = 1, this%nb_total_nodes
                    if(this%nodes(i_node)%keep) then
                        print*, "Node: ", i_node, this%nodes(i_node)%id, this%nodes(i_node)%coor
                    end if
                end do
            end if
            ASSERT(this%nb_nodes <= this%nb_total_nodes)
    end subroutine
!
! ==================================================================================================
!
    subroutine check_mesh(this)
!
        implicit none
!
            class(Mmesh), intent(in) :: this
! -----------------------------------------------------------------------
            integer :: i_cell, i_node, i_edge, i_face, i_volume, nno, nno1, nno2
            integer :: nb_edges, edges_type(12), edges_loc(3,12), edge_id
            integer :: face_type, volu_type, face_id
            integer :: nb_faces, faces_type(6), faces_loc(9,6)
!
! --- Check Nodes
            do i_cell = 1, this%nb_cells
                nno = this%converter%nno(this%cells(i_cell)%type)
                do i_node = 1, nno
                    ASSERT(this%nodes(this%cells(i_cell)%nodes(i_node))%keep)
                end do
            end do
!
! --- Check Edges
            do i_edge = 1, this%nb_edges
                nno = this%converter%nno(this%edges(i_edge)%type)
                do i_node = 1, nno
                    ASSERT(this%edges(i_edge)%nodes(i_node) > 0)
                end do
            end do
!
! --- Check Faces
            do i_face = 1, this%nb_faces
                face_type = this%faces(i_face)%type
                nno = this%converter%nno(face_type)
                do i_node = 1, nno
                    ASSERT(this%faces(i_face)%nodes(i_node) > 0)
                end do
!
                call numbering_edge(face_type, nb_edges, edges_type, edges_loc)
                ASSERT(nb_edges == this%faces(i_face)%nb_edges)
                do i_edge = 1, this%faces(i_face)%nb_edges
                    edge_id = this%faces(i_face)%edges(i_edge)
                    nno1 = this%converter%nno(this%edges(edge_id)%type)
                    nno2 = this%converter%nno(edges_type(i_edge))
                    ASSERT(nno1 >= nno2)
                end do
            end do
!
! --- Check Volumes
            do i_volume = 1, this%nb_volumes
                volu_type = this%volumes(i_volume)%type
                nno = this%converter%nno(volu_type)
                do i_node = 1, nno
                    ASSERT(this%volumes(i_volume)%nodes(i_node) > 0)
                end do
!
                call numbering_edge(volu_type, nb_edges, edges_type, edges_loc)
                ASSERT(nb_edges == this%volumes(i_volume)%nb_edges)
                do i_edge = 1, this%volumes(i_volume)%nb_edges
                    edge_id = this%volumes(i_volume)%edges(i_edge)
                    nno1 = this%converter%nno(this%edges(edge_id)%type)
                    nno2 = this%converter%nno(edges_type(i_edge))
                    ASSERT(nno1 >= nno2)
                end do
!
                call numbering_face(volu_type, nb_faces, faces_type, faces_loc)
                ASSERT(nb_faces == this%volumes(i_volume)%nb_faces)
                do i_face = 1, this%volumes(i_volume)%nb_faces
                    face_id = this%volumes(i_volume)%faces(i_face)
                    nno1 = this%converter%nno(this%faces(face_id)%type)
                    nno2 = this%converter%nno(faces_type(i_face))
                    ASSERT(nno1 >= nno2)
                end do
            end do
    end subroutine
!
! ==================================================================================================
!
    subroutine copy_group_no(this, grpnoe, gpptnn)
!
        implicit none
!
            class(Mmesh), intent(in) :: this
            character(len=24), intent(in) :: grpnoe, gpptnn
! -----------------------------------------------------------------------
            integer :: i_node, nb_nodes_in, nb_nodes_out, codret, nb_grno_out
            integer :: i_group, node_id, nb_grno_in
            character(len=24) :: grno_in, nomgrp
            integer, pointer :: nodes_in(:) => null()
            integer, pointer :: grno_out(:) => null()
            integer, pointer :: nodes_out(:) => null()
!
            call jemarq()
!
            grno_in = this%mesh_in//'.GROUPENO'
!
            call jedetr(grpnoe)
            call jedetr(gpptnn)
            call jeexin(grno_in, codret)
            if (codret .eq. 0) goto 999

            call jelira(grno_in, 'NOMUTI', nb_grno_in)
            AS_ALLOCATE(vi=grno_out, size=nb_grno_in)
            grno_out(:) = 0
            nb_grno_out = 0
!
! --- Find groups
            do i_group = 1, nb_grno_in
                call jeveuo(jexnum(grno_in, i_group), 'L', vi=nodes_in)
                call jelira(jexnum(grno_in, i_group), 'LONUTI', nb_nodes_in)
                do i_node = 1, nb_nodes_in
                    node_id = nodes_in(i_node)
                    if(this%nodes(node_id)%keep) then
                        grno_out(i_group) = grno_out(i_group) + 1
                    end if
                end do
                if(grno_out(i_group) > 0) then
                    nb_grno_out = nb_grno_out + 1
                end if
            end do
!
! --- Create groups
            if(nb_grno_out > 0) then
                call jecreo(gpptnn, 'G N K24')
                call jeecra(gpptnn, 'NOMMAX', nb_grno_out)
                call jecrec(grpnoe, 'G V I', 'NO '//gpptnn, 'DISPERSE', 'VARIABLE', nb_grno_out)
!
                do i_group = 1, nb_grno_in
                    if(grno_out(i_group) > 0) then
                        call jenuno(jexnum(grno_in, i_group), nomgrp)
                        call jecroc(jexnom(grpnoe, nomgrp))
                        call jeveuo(jexnum(grno_in, i_group), 'L', vi=nodes_in)
                        call jelira(jexnum(grno_in, i_group), 'LONUTI', nb_nodes_in)
                        call jeecra(jexnom(grpnoe, nomgrp), 'LONMAX', grno_out(i_group))
                        call jeecra(jexnom(grpnoe, nomgrp), 'LONUTI', grno_out(i_group))
                        call jeveuo(jexnom(grpnoe, nomgrp), 'E', vi=nodes_out)
                        nb_nodes_out = 0
                        do i_node = 1, nb_nodes_in
                            node_id = nodes_in(i_node)
                            if(this%nodes(node_id)%keep) then
                                nb_nodes_out = nb_nodes_out + 1
                                nodes_out(nb_nodes_out) = this%nodes(node_id)%id
                            end if
                        end do
                    end if
                end do
            end if
!
            AS_DEALLOCATE(vi=grno_out)
!
999 continue
!
            call jedema()
!
    end subroutine
!
    !
! ==================================================================================================
!
    subroutine increase_memory(this, object, new_size)
!
        implicit none
!
            class(Mmesh), intent(inout) :: this
            character(len=8), intent(in) :: object
            integer, intent(in) :: new_size
! -----------------------------------------------------------------------
            integer :: old_size
            type(Mnode), allocatable :: nodes(:)
            type(Medge), allocatable :: edges(:)
            type(Mface), allocatable :: faces(:)
            type(Mvolume), allocatable :: volumes(:)
            type(Mcell), allocatable :: cells(:)
!
            if(object == "NODES") then
                old_size = size(this%nodes)
                allocate(nodes(old_size))
                nodes(1:old_size) = this%nodes(1:old_size)
                deallocate(this%nodes)
                allocate(this%nodes(new_size))
                this%nodes(1:old_size) = nodes(1:old_size)
                deallocate(nodes)
                this%max_nodes = new_size
            elseif(object == "EDGES") then
                old_size = size(this%edges)
                allocate(edges(old_size))
                edges(1:old_size) = this%edges(1:old_size)
                deallocate(this%edges)
                allocate(this%edges(new_size))
                this%edges(1:old_size) = edges(1:old_size)
                deallocate(edges)
                this%max_edges = new_size
            elseif(object == "FACES") then
                old_size = size(this%faces)
                allocate(faces(old_size))
                faces(1:old_size) = this%faces(1:old_size)
                deallocate(this%faces)
                allocate(this%faces(new_size))
                this%faces(1:old_size) = faces(1:old_size)
                deallocate(faces)
                this%max_faces = new_size
            elseif(object == "VOLUMES") then
                old_size = size(this%volumes)
                allocate(volumes(old_size))
                volumes(1:old_size) = this%volumes(1:old_size)
                deallocate(this%volumes)
                allocate(this%volumes(new_size))
                this%volumes(1:old_size) = volumes(1:old_size)
                deallocate(volumes)
                this%max_volumes = new_size
            elseif(object == "CELLS") then
                old_size = size(this%cells)
                allocate(cells(old_size))
                cells(1:old_size) = this%cells(1:old_size)
                deallocate(this%cells)
                allocate(this%cells(new_size))
                this%cells(1:old_size) = cells(1:old_size)
                deallocate(cells)
                this%max_cells = new_size
            else
                ASSERT(ASTER_FALSE)
            end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine create_joints(this, mesh_out)
!
        implicit none
!
            class(Mmesh), intent(in) :: this
            character(len=8), intent(in) :: mesh_out
! ------------------------------------------------------------------
!
        if(isParallelMesh(mesh_out)) then
            ASSERT(ASTER_FALSE)
! --- TO DO
!
            ASSERT(this%nb_nodes > 0)
! --- 1: Create new global numbering .NULOGL
!
! --- 2: To known who is the own the nodes .NOEX
!
! --- 3: Create joint between subdomains .DOMJOINT
        end if
    end subroutine
!
end module
