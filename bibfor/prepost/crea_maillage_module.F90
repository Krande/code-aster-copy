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
! person_in_charge: nicolas.pignet at edf.fr
!
module crea_maillage_module
!
    use sort_module
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_sendrecv_i.h"
#include "asterc/asmpi_sendrecv_r.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/build_tree_comm.h"
#include "asterfort/codent.h"
#include "asterfort/elraga.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/isParallelMesh.h"
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
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/sdmail.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! Module to create or modify a mesh.
!
! --------------------------------------------------------------------------------------------------
!
!
    type Mconverter
        aster_logical :: to_convert(MT_NTYMAX) = ASTER_FALSE
        integer(kind=8) :: convert_to(MT_NTYMAX) = 0
        integer(kind=8) :: convert_max(MT_NTYMAX) = 0
        character(len=8) :: name(MT_NTYMAX) = ' '
        character(len=8) :: short_name(MT_NTYMAX) = ' '
        integer(kind=8) :: cata_type(MT_NTYMAX) = 0
        integer(kind=8) :: map_type(MT_NTYMAX) = 0
        integer(kind=8) :: dim(MT_NTYMAX) = 0
        integer(kind=8) :: nno(MT_NTYMAX) = 0
        integer(kind=8) :: nnos(MT_NTYMAX) = 0
! ----- member functions
    contains
        procedure, public, pass :: init => init_conv
        procedure, public, pass :: add_conversion
    end type
!
    type Mcell
        integer(kind=8) :: type = 0
        integer(kind=8) :: dim = 0
        integer(kind=8) :: id = 0
        integer(kind=8) :: ss_id = 0
        integer(kind=8) :: nodes(27) = 0
        integer(kind=8) :: child(10) = 0
        integer(kind=8) :: nb_child = 0
        aster_logical :: keep = ASTER_TRUE
    end type
!
    type Mvolume
        integer(kind=8) :: type = 0
        integer(kind=8) :: nodes(27) = 0
        integer(kind=8) :: nb_faces = 0, faces(6) = 0
        integer(kind=8) :: nb_edges = 0, edges(12) = 0
        integer(kind=8) :: parent = -1
        integer(kind=8) :: isub = 0
        integer(kind=8) :: cell_id = 0
    end type
!
    type Mface
        integer(kind=8) :: type = 0
        integer(kind=8) :: nnos = 0
        integer(kind=8) :: nodes(9) = 0
        integer(kind=8) :: nno_sort(9) = 0
        integer(kind=8) :: nb_volumes = 0, volumes(2) = 0
        integer(kind=8) :: nb_edges = 0, edges(4) = 0
        integer(kind=8) :: parent = -1
        integer(kind=8) :: isub = 0
        integer(kind=8) :: cell_id = 0
    end type
!
    type Medge
        integer(kind=8) :: type = 0
        integer(kind=8) :: nodes(4) = 0
        integer(kind=8) :: nno_sort(4) = 0
        integer(kind=8) :: nb_faces = 0, faces(20) = 0
        integer(kind=8) :: nb_volumes = 0, volumes(25) = 0
        integer(kind=8) :: parent = -1
        integer(kind=8) :: isub = 0
        integer(kind=8) :: cell_id = 0
    end type
!
    type Mnode
        integer(kind=8) :: id
        real(kind=8) :: coor(3) = 0.d0
        aster_logical :: keep = ASTER_FALSE
        aster_logical :: orphelan = ASTER_FALSE
        integer(kind=8) :: owner = -1
! used to improve search of edges and faces
! it could be improved a lot to decrease memory consumption
        integer(kind=8) :: max_faces = 0, nb_faces = 0
        integer(kind=8) :: max_edges = 0, nb_edges = 0
        integer(kind=8), allocatable :: faces(:)
        integer(kind=8), allocatable :: edges(:)
    end type
!
    type Mmesh
        integer(kind=8) :: nb_nodes = 0, nb_edges = 0, nb_faces = 0, nb_volumes = 0, nb_cells = 0
        integer(kind=8) :: nb_total_nodes = 0, nb_total_cells = 0
        integer(kind=8) :: max_nodes = 0, max_edges = 0, max_faces = 0, max_volumes = 0
        integer(kind=8) :: max_cells = 0, dim_mesh = 0, nb_layer = 0, lastLayerSize = 0
        integer(kind=8) :: nb_edges_dege = 0, nb_faces_dege = 0
        integer(kind=8) :: max_faces_dege = 0, max_edges_dege = 0
        integer(kind=8), allocatable :: faces_dege(:), edges_dege(:)
!
        integer(kind=4), allocatable :: lastGhostsLayer(:)

        type(Mnode), allocatable :: nodes(:)
        type(Medge), allocatable :: edges(:)
        type(Mface), allocatable :: faces(:)
        type(Mvolume), allocatable :: volumes(:)
        type(Mcell), allocatable :: cells(:)
        type(Mconverter) :: converter
!
        integer(kind=8) :: nb_level = 0
!
        character(len=8) :: mesh_in = ' '
        integer(kind=8), pointer :: connex_in(:) => null()
        integer(kind=8), pointer :: connex_map_in(:) => null()
!
        integer(kind=8), pointer :: v_typema(:) => null()
        aster_logical :: debug = ASTER_FALSE
        aster_logical :: isHPC
        integer(kind=8) :: info = 0
! ----- member functions
    contains
        procedure, public, pass :: add_initial_cell
        procedure, public, pass :: check_mesh
        procedure, public, pass :: clean => clean_mesh
        procedure, public, pass :: convert_cells
        procedure, public, pass :: copy_mesh
        procedure, public, pass :: create_joints
        procedure, public, pass :: init => init_mesh
        procedure, public, pass :: fix => fix_mesh
        procedure, public, pass :: refine
        procedure, private, pass :: add_edge
        procedure, private, pass :: add_face
        procedure, private, pass :: add_node
        procedure, private, pass :: add_point1
        procedure, private, pass :: add_volume
        procedure, private, pass :: barycenter
        procedure, private, pass :: barycenter_edge
        procedure, private, pass :: barycenter_face
        procedure, private, pass :: barycenter_volume
        procedure, private, pass :: barycenter_fe
        procedure, private, pass :: convert_edge
        procedure, private, pass :: convert_face
        procedure, private, pass :: convert_volume
        procedure, private, pass :: copy_group_ma
        procedure, private, pass :: copy_group_no
        procedure, private, pass :: find_edge
        procedure, private, pass :: find_face
        procedure, private, pass :: increase_memory
        procedure, private, pass :: numbering_nodes
        procedure, private, pass :: owner_cell
        procedure, private, pass :: refine_cell
        procedure, private, pass :: sub_cells
        procedure, private, pass :: update
    end type
!
!===================================================================================================
!
!===================================================================================================
!
    public :: Medge, Mface, Mcell, Mmesh, Mconverter
    private :: numbering_edge, numbering_face, dividing_cell, mult_elem
    private :: sort_nodes_face, fix_measure
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
        integer(kind=8) :: from_i, to_i
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
    subroutine find_edge(this, nno, nno_sort, find, edge_id)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: nno_sort(3), nno
        aster_logical, intent(out) :: find
        integer(kind=8), intent(out) :: edge_id
!
! Find edge_id of edges
!
        integer(kind=8) :: i_edge, nb_edges, edge_i, edge_i_error
        aster_logical :: ok
!
        find = ASTER_FALSE
        edge_id = 0
        edge_i_error = 0
!
        nb_edges = this%nodes(nno_sort(1))%nb_edges
        if (nb_edges > 0) then
            do i_edge = 1, nb_edges
                edge_i = this%nodes(nno_sort(1))%edges(i_edge)
                ok = ASTER_TRUE
                if (nno_sort(1) .ne. this%edges(edge_i)%nno_sort(1)) then
                    ok = ASTER_FALSE
                else
                    if (nno_sort(2) .ne. this%edges(edge_i)%nno_sort(2)) then
                        ok = ASTER_FALSE
                    end if
                end if
                if (ok) then
                    if (nno_sort(nno) .ne. this%edges(edge_i)%nno_sort(nno)) then
                        edge_i_error = edge_i
                        ok = ASTER_FALSE
                    end if
                end if
                if (ok) then
                    find = ASTER_TRUE
                    edge_id = edge_i
                    exit
                end if
            end do
        end if
!
        if (.not. find .and. edge_i_error > 0) then
            this%nb_edges_dege = this%nb_edges_dege+1
            if (2*this%nb_edges_dege > this%max_edges_dege) then
                call this%increase_memory("EDGES_DEG", 2*this%max_edges_dege)
            end if
            this%edges_dege(2*(this%nb_edges_dege-1)+1) = edge_i_error
            this%edges_dege(2*(this%nb_edges_dege-1)+2) = this%nb_edges+1
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
    subroutine find_face(this, nno, nnos, nno_sort, find, face_id)
!
        implicit none
!
        class(Mmesh), intent(inout) ::this
        integer(kind=8), intent(in) :: nnos, nno
        integer(kind=8), intent(in) :: nno_sort(9)
        aster_logical, intent(out) :: find
        integer(kind=8), intent(out) :: face_id
!
! Find face_id of face
!
        integer(kind=8) :: i_face, i_node, nb_faces, face_i, face_i_error
        aster_logical :: ok
!
        find = ASTER_FALSE
        face_id = 0
        face_i_error = 0
!
        nb_faces = this%nodes(nno_sort(1))%nb_faces
        if (nb_faces > 0) then
            do i_face = 1, nb_faces
                face_i = this%nodes(nno_sort(1))%faces(i_face)
                if (nnos == this%faces(face_i)%nnos) then
                    ok = ASTER_TRUE
                    do i_node = 1, nnos
                        if (nno_sort(i_node) .ne. this%faces(face_i)%nno_sort(i_node)) then
                            ok = ASTER_FALSE
                            exit
                        end if
                    end do
                else
                    ok = ASTER_FALSE
                end if
                if (ok) then
                    do i_node = nnos+1, nno
                        if (nno_sort(i_node) .ne. this%faces(face_i)%nno_sort(i_node)) then
                            face_i_error = face_i
                            ok = ASTER_FALSE
                            exit
                        end if
                    end do
                end if
                if (ok) then
                    find = ASTER_TRUE
                    face_id = face_i
                    exit
                end if
            end do
        end if
!
        if (.not. find .and. face_i_error > 0) then
            this%nb_faces_dege = this%nb_faces_dege+1
            if (2*this%nb_faces_dege > this%max_faces_dege) then
                call this%increase_memory("FACES_DEG", 2*this%max_faces_dege)
            end if
            this%faces_dege(2*(this%nb_faces_dege-1)+1) = face_i_error
            this%faces_dege(2*(this%nb_faces_dege-1)+2) = this%nb_faces+1
        end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine sort_nodes_face(nno, nnos, nodes)
!
        implicit none
!
        integer(kind=8), intent(in) :: nno, nnos
        integer(kind=8), intent(inout) :: nodes(1:nno)
!
        ASSERT(nno <= 9)
!
        call qsort(nodes(1:nnos))
!
        if (nno > nnos) then
            call qsort(nodes(nnos+1:2*nnos))
        end if
!
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
        integer(kind=8) :: i_type, type_nume
!
        call jemarq()
!
        this%name(MT_POI1) = "POI1"
        this%name(MT_SEG2) = "SEG2"
        this%name(MT_SEG3) = "SEG3"
        this%name(MT_SEG4) = "SEG4"
        this%name(MT_TRIA3) = "TRIA3"
        this%name(MT_TRIA6) = "TRIA6"
        this%name(MT_TRIA7) = "TRIA7"
        this%name(MT_QUAD4) = "QUAD4"
        this%name(MT_QUAD8) = "QUAD8"
        this%name(MT_QUAD9) = "QUAD9"
        this%name(MT_TETRA4) = "TETRA4"
        this%name(MT_TETRA10) = "TETRA10"
        this%name(MT_TETRA15) = "TETRA15"
        this%name(MT_HEXA8) = "HEXA8"
        this%name(MT_HEXA9) = "HEXA9"
        this%name(MT_HEXA20) = "HEXA20"
        this%name(MT_HEXA27) = "HEXA27"
        this%name(MT_PYRAM5) = "PYRAM5"
        this%name(MT_PYRAM13) = "PYRAM13"
        this%name(MT_PYRAM19) = "PYRAM19"
        this%name(MT_PENTA6) = "PENTA6"
        this%name(MT_PENTA7) = "PENTA7"
        this%name(MT_PENTA15) = "PENTA15"
        this%name(MT_PENTA18) = "PENTA18"
        this%name(MT_PENTA21) = "PENTA21"
!
        this%short_name(MT_POI1) = "PO1"
        this%short_name(MT_SEG2) = "SE2"
        this%short_name(MT_SEG3) = "SE3"
        this%short_name(MT_SEG4) = "SE4"
        this%short_name(MT_TRIA3) = "TR3"
        this%short_name(MT_TRIA6) = "TR6"
        this%short_name(MT_TRIA7) = "TR7"
        this%short_name(MT_QUAD4) = "QU4"
        this%short_name(MT_QUAD8) = "QU8"
        this%short_name(MT_QUAD9) = "QU9"
        this%short_name(MT_TETRA4) = "TE4"
        this%short_name(MT_TETRA10) = "T10"
        this%short_name(MT_TETRA15) = "T15"
        this%short_name(MT_HEXA8) = "HE8"
        this%short_name(MT_HEXA9) = "HE9"
        this%short_name(MT_HEXA20) = "H20"
        this%short_name(MT_HEXA27) = "H27"
        this%short_name(MT_PYRAM5) = "PY5"
        this%short_name(MT_PYRAM13) = "P13"
        this%short_name(MT_PYRAM19) = "P19"
        this%short_name(MT_PENTA6) = "PE6"
        this%short_name(MT_PENTA7) = "PE7"
        this%short_name(MT_PENTA15) = "P15"
        this%short_name(MT_PENTA18) = "P18"
        this%short_name(MT_PENTA21) = "P21"
!
        this%convert_max(MT_POI1) = MT_POI1
        this%convert_max(MT_SEG2) = MT_SEG3
        this%convert_max(MT_SEG3) = MT_SEG3
        this%convert_max(MT_SEG4) = MT_SEG4
        this%convert_max(MT_TRIA3) = MT_TRIA7
        this%convert_max(MT_TRIA6) = MT_TRIA7
        this%convert_max(MT_TRIA7) = MT_TRIA7
        this%convert_max(MT_QUAD4) = MT_QUAD9
        this%convert_max(MT_QUAD8) = MT_QUAD9
        this%convert_max(MT_QUAD9) = MT_QUAD9
        this%convert_max(MT_TETRA4) = MT_TETRA15
        this%convert_max(MT_TETRA10) = MT_TETRA15
        this%convert_max(MT_TETRA15) = MT_TETRA15
        this%convert_max(MT_HEXA8) = MT_HEXA27
        this%convert_max(MT_HEXA9) = MT_HEXA27
        this%convert_max(MT_HEXA20) = MT_HEXA27
        this%convert_max(MT_HEXA27) = MT_HEXA27
        this%convert_max(MT_PYRAM5) = MT_PYRAM19
        this%convert_max(MT_PYRAM13) = MT_PYRAM19
        this%convert_max(MT_PYRAM19) = MT_PYRAM19
        this%convert_max(MT_PENTA6) = MT_PENTA21
        this%convert_max(MT_PENTA7) = MT_PENTA21
        this%convert_max(MT_PENTA15) = MT_PENTA21
        this%convert_max(MT_PENTA18) = MT_PENTA21
        this%convert_max(MT_PENTA21) = MT_PENTA21
!
        do i_type = 1, MT_NTYMAX
            if (this%name(i_type) .ne. ' ') then
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
        integer(kind=8), optional :: info
! --------------------------------------------------------------------------------------------------
! The idea is to read a given mesh (mesh_in). Internally, all the cells are stored with all nodes
! possibles for a cell. The internal cells stored are POI1, SEG3, SEG4, TRIA7, QUAD9,
! TETRA15, HEXA27, PENTA21 and PYRAM19
! For an other cell type is only necessary to know which nodes to use
! --------------------------------------------------------------------------------------------------
        character(len=19) :: joints
        integer(kind=8), pointer :: v_mesh_dime(:) => null()
        integer(kind=8), pointer :: v_noex(:) => null()
        integer(kind=8), pointer :: v_nblg(:) => null()
        real(kind=8), pointer :: v_coor(:) => null()
        integer(kind=8) :: nb_elem_mesh, nb_node_mesh, i_node, iret
        integer(kind=8) :: i_cell, nno, node_id, owner, i_type, cell_type
        real(kind=8):: start, end
        integer(kind=8), allocatable, dimension(:, :) :: list_type_cells
        integer(kind=8), parameter :: nb_type = 10
        integer(kind=8) :: nb_type_cells(nb_type)
!
        call jemarq()
        call this%converter%init()
!
        if (present(info)) then
            this%info = info
        end if
!
        if (this%info >= 2) then
            print *, "Creating mesh..."
            call cpu_time(start)
        end if
!
        this%isHPC = isParallelMesh(mesh_in)
!
        call jeveuo(mesh_in//'.DIME', 'L', vi=v_mesh_dime)
        this%dim_mesh = v_mesh_dime(6)
        nb_elem_mesh = v_mesh_dime(3)
        nb_node_mesh = v_mesh_dime(1)
        this%mesh_in = mesh_in
!
        call jeveuo(mesh_in//'.CONNEX', 'L', vi=this%connex_in)
        call jeveuo(jexatr(mesh_in//'.CONNEX', 'LONCUM'), 'L', vi=this%connex_map_in)
!
        this%max_cells = max(1, nb_elem_mesh)
        this%max_volumes = max(1, nb_elem_mesh)
        this%max_faces = max(1, 6*nb_elem_mesh)
        this%max_edges = max(1, 12*nb_elem_mesh)
        this%max_nodes = max(1, 5*nb_node_mesh)
        this%max_faces_dege = 1000
        this%max_edges_dege = 1000
!
        allocate (this%cells(this%max_cells))
        allocate (this%volumes(this%max_volumes))
        allocate (this%faces(this%max_faces))
        allocate (this%edges(this%max_edges))
        allocate (this%nodes(this%max_nodes))
        allocate (this%lastGhostsLayer(this%max_nodes))
        allocate (this%faces_dege(this%max_faces_dege))
        allocate (this%edges_dege(this%max_edges_dege))
!
        this%faces_dege = 0
        this%edges_dege = 0
!
        call jeveuo(mesh_in//'.TYPMAIL', 'L', vi=this%v_typema)
        call jeveuo(mesh_in//'.COORDO    .VALE', 'L', vr=v_coor)
        if (this%isHPC) then
            call jeveuo(mesh_in//'.NOEX', 'L', vi=v_noex)
            joints = mesh_in//".JOIN"
            call jeexin(joints//'.NBLG', iret)
            if (iret > 0) then
                call jeveuo(joints//'.NBLG', 'L', vi=v_nblg)
                this%nb_layer = v_nblg(1)
                ASSERT(this%nb_layer > 0)
            end if
        end if
!
! --- Fill mesh
!
        owner = 0
        do i_node = 1, nb_node_mesh
            if (this%isHPC) then
                owner = v_noex(i_node)
            end if
            node_id = this%add_node(v_coor(3*(i_node-1)+1:3*(i_node-1)+3), owner)
            ASSERT(i_node == node_id)
            this%nodes(node_id)%orphelan = ASTER_TRUE
            this%nodes(node_id)%max_faces = 30
            allocate (this%nodes(node_id)%faces(this%nodes(node_id)%max_faces))
            this%nodes(node_id)%max_edges = 30
            allocate (this%nodes(node_id)%edges(this%nodes(node_id)%max_edges))
        end do
!
! --- Read cells
        allocate (list_type_cells(nb_type, nb_elem_mesh))
!
!
! --- Order to read cell - hard coded rule
        nb_type_cells = 0
        do i_cell = 1, nb_elem_mesh
            cell_type = this%converter%map_type(this%v_typema(i_cell))
!
            select case (cell_type)
            case (MT_POI1)
                nb_type_cells(1) = nb_type_cells(1)+1
                list_type_cells(1, nb_type_cells(1)) = i_cell
            case (MT_SEG4)
                nb_type_cells(2) = nb_type_cells(2)+1
                list_type_cells(2, nb_type_cells(2)) = i_cell
            case (MT_SEG3)
                nb_type_cells(3) = nb_type_cells(3)+1
                list_type_cells(3, nb_type_cells(3)) = i_cell
            case (MT_TRIA7, MT_QUAD9)
                nb_type_cells(4) = nb_type_cells(4)+1
                list_type_cells(4, nb_type_cells(4)) = i_cell
            case (MT_HEXA27, MT_PENTA21, MT_PYRAM19, MT_TETRA15)
                nb_type_cells(5) = nb_type_cells(5)+1
                list_type_cells(5, nb_type_cells(5)) = i_cell
            case (MT_TRIA6, MT_QUAD8)
                nb_type_cells(6) = nb_type_cells(6)+1
                list_type_cells(6, nb_type_cells(6)) = i_cell
            case (MT_PENTA18, MT_HEXA20, MT_PENTA15, MT_PYRAM13, MT_TETRA10)
                nb_type_cells(7) = nb_type_cells(7)+1
                list_type_cells(7, nb_type_cells(7)) = i_cell
            case (MT_SEG2)
                nb_type_cells(8) = nb_type_cells(8)+1
                list_type_cells(8, nb_type_cells(8)) = i_cell
            case (MT_TRIA3, MT_QUAD4)
                nb_type_cells(9) = nb_type_cells(9)+1
                list_type_cells(9, nb_type_cells(9)) = i_cell
            case (MT_HEXA9, MT_PENTA7, MT_HEXA8, MT_PENTA6, MT_PYRAM5, MT_TETRA4)
                nb_type_cells(10) = nb_type_cells(10)+1
                list_type_cells(10, nb_type_cells(10)) = i_cell
            case default
                ASSERT(ASTER_FALSE)
            end select
        end do
!
        do i_type = 1, nb_type
            do i_cell = 1, nb_type_cells(i_type)
                call this%add_initial_cell(list_type_cells(i_type, i_cell))
            end do
        end do

        deallocate (list_type_cells)
!
! --- Search orphelan nodes - to keep at the end
        this%nb_total_cells = this%nb_cells
        do i_cell = 1, this%nb_cells
            nno = this%converter%nno(this%cells(i_cell)%type)
            do i_node = 1, nno
                node_id = this%cells(i_cell)%nodes(i_node)
                this%nodes(node_id)%orphelan = ASTER_FALSE
            end do
        end do
!
        if (this%info >= 2) then
            call cpu_time(end)
            print *, "... in ", end-start, " seconds."
        end if
!
        call jedema()
!
    end subroutine
!
! ==================================================================================================
!
    subroutine fix_mesh(this, remove_orphelan, outward_normal, positive_measure)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        aster_logical, intent(in) :: remove_orphelan, outward_normal, positive_measure
!
        integer(kind=8) :: i_edge, e1, e2, n1, n2, n3, i_face, f1, f2
        integer(kind=8) :: i, j, i_node, nno, i_volu, v_id, nnos, c_id
        integer(kind=8), allocatable :: map_nodes(:)
        real(kind=8) :: new_coor(3), end, start, v0(3), v1(3), no(3)
        real(kind=8) :: bar_v(3), bar_f(3), nvf(3)

!
        if (this%info >= 2) then
            print *, "Fixing mesh..."
            call cpu_time(start)
        end if
!
        allocate (map_nodes(this%nb_total_nodes))
        do i_node = 1, this%nb_total_nodes
            map_nodes(i_node) = i_node
        end do
!
        do i_edge = 1, this%nb_edges_dege
            e1 = this%edges_dege(2*(i_edge-1)+1)
            e2 = this%edges_dege(2*(i_edge-1)+2)
            n1 = this%edges(e1)%nno_sort(3)
            n2 = this%edges(e2)%nno_sort(3)
            ASSERT(n1*n2 > 0)
!
!           use barycenter of two middle nodes
            new_coor = 0.5d0*(this%nodes(n1)%coor+this%nodes(n2)%coor)
            this%nodes(n1)%coor = new_coor
            this%nodes(n2)%keep = ASTER_FALSE
            map_nodes(n2) = n1
            this%nb_nodes = this%nb_nodes-1
            if (this%edges(e2)%cell_id > 0) then
                this%cells(this%edges(e2)%cell_id)%keep = ASTER_FALSE
            end if
!
            do i_face = 1, this%edges(e2)%nb_faces
                f2 = this%edges(e2)%faces(i_face)
                do i = 1, this%faces(f2)%nb_edges
                    if (e2 == this%faces(f2)%edges(i)) then
                        this%faces(f2)%edges(i) = e1
                        nno = this%converter%nno(this%faces(f2)%type)
                        do i_node = 1, nno
                            if (n2 == this%faces(f2)%nodes(i_node)) then
                                this%faces(f2)%nodes(i_node) = n1
                                do j = 1, nno
                                    if (n2 == this%faces(f2)%nno_sort(j)) then
                                        this%faces(f2)%nno_sort(j) = n1
                                        exit
                                    end if
                                end do
                                exit
                            end if
                        end do
                        exit
                    end if
                end do
            end do
!
            do i_volu = 1, this%edges(e2)%nb_volumes
                v_id = this%edges(e2)%volumes(i_volu)
                do i = 1, this%volumes(v_id)%nb_edges
                    if (e2 == this%volumes(v_id)%edges(i)) then
                        this%volumes(v_id)%edges(i) = e1
                        nno = this%converter%nno(this%volumes(v_id)%type)
                        do i_node = 1, nno
                            if (n2 == this%volumes(v_id)%nodes(i_node)) then
                                this%volumes(v_id)%nodes(i_node) = n1
                                exit
                            end if
                        end do
                        exit
                    end if
                    exit
                end do
            end do
!
        end do
!
        do i_face = 1, this%nb_faces_dege
            f1 = this%faces_dege(2*(i_face-1)+1)
            f2 = this%faces_dege(2*(i_face-1)+2)
            nno = this%converter%nno(this%faces(f1)%type)

            n1 = this%faces(f1)%nodes(nno)
            n2 = this%faces(f2)%nodes(nno)
            ASSERT(n1*n2 > 0)
!
!           use barycenter of two middle nodes
            new_coor = 0.5d0*(this%nodes(n1)%coor+this%nodes(n2)%coor)
            this%nodes(n1)%coor = new_coor
            this%nodes(n2)%keep = ASTER_FALSE
            map_nodes(n2) = n1
            this%nb_nodes = this%nb_nodes-1
            if (this%faces(f2)%cell_id > 0) then
                this%cells(this%faces(f2)%cell_id)%keep = ASTER_FALSE
            end if
            ASSERT(this%faces(f1)%nb_volumes == 1)
            ASSERT(this%faces(f2)%nb_volumes == 1)
            this%faces(f1)%nb_volumes = 2
            this%faces(f1)%volumes(2) = this%faces(f2)%volumes(1)
!
            do i_volu = 1, this%faces(f2)%nb_volumes
                v_id = this%faces(f2)%volumes(i_volu)
                do i = 1, this%volumes(v_id)%nb_faces
                    if (f2 == this%volumes(v_id)%faces(i)) then
                        this%volumes(v_id)%faces(i) = f1
                        nno = this%converter%nno(this%volumes(v_id)%type)
                        do i_node = 1, nno
                            if (n2 == this%volumes(v_id)%nodes(i_node)) then
                                this%volumes(v_id)%nodes(i_node) = n1
                                exit
                            end if
                        end do
                        exit
                    end if
                    exit
                end do
            end do
!
        end do
!
! --- Renumbering nodes
        do i = 1, this%nb_cells
            if (this%cells(i)%keep) then
                nno = this%converter%nno(this%cells(i)%type)
                do i_node = 1, nno
                    n2 = this%cells(i)%nodes(i_node)
                    n1 = map_nodes(n2)
                    this%cells(i)%nodes(i_node) = n1
                end do
            end if
        end do
!
        deallocate (map_nodes)
!
        if (remove_orphelan) then
            do i_node = 1, this%nb_nodes
                if (this%nodes(i_node)%orphelan) then
                    this%nodes(i_node)%keep = ASTER_FALSE
                    this%nodes(i_node)%orphelan = ASTER_FALSE
                end if
            end do
        end if
!
        if (outward_normal) then
            do i_face = 1, this%nb_faces
                if (this%faces(i_face)%nb_volumes == 1) then
                    nnos = this%converter%nnos(this%faces(i_face)%type)
                    n1 = this%faces(i_face)%nodes(1)
                    n2 = this%faces(i_face)%nodes(2)
                    n3 = this%faces(i_face)%nodes(nnos)
                    v0 = this%nodes(n2)%coor-this%nodes(n1)%coor
                    v1 = this%nodes(n3)%coor-this%nodes(n1)%coor
                    no(1) = v0(2)*v1(3)-v0(3)*v1(2)
                    no(2) = v0(3)*v1(1)-v0(1)*v1(3)
                    no(3) = v0(1)*v1(2)-v0(2)*v1(1)
                    no = no/norm2(no)
                    bar_f = this%barycenter_face(i_face)
                    bar_v = this%barycenter_volume(this%faces(i_face)%volumes(1))
                    nvf = bar_f-bar_v
                    nvf = nvf/norm2(nvf)
                    if (dot_product(nvf, no) <= 0.d0) then
                        this%faces(i_face)%nodes(2) = n3
                        this%faces(i_face)%nodes(nnos) = n2
                        n2 = this%faces(i_face)%nodes(nnos+1)
                        n3 = this%faces(i_face)%nodes(2*nnos)
                        this%faces(i_face)%nodes(nnos+1) = n3
                        this%faces(i_face)%nodes(2*nnos) = n2
                        if (nnos == 4) then
                            n2 = this%faces(i_face)%nodes(nnos+2)
                            n3 = this%faces(i_face)%nodes(nnos+3)
                            this%faces(i_face)%nodes(nnos+2) = n3
                            this%faces(i_face)%nodes(nnos+3) = n2
                        end if
                        c_id = this%faces(i_face)%cell_id
                        if (c_id > 0) then
                            nno = this%converter%nno(this%cells(c_id)%type)
                            this%cells(c_id)%nodes(1:nno) = this%faces(i_face)%nodes(1:nno)
                        end if
                    end if
                end if
            end do
        end if
!
        if (positive_measure) then
            do i_volu = 1, this%nb_volumes
                call fix_measure(this, i_volu, 3)
            end do
        end if
!
! --- Keep only necessary nodes
        call this%update()
!
        if (this%info >= 2) then
            call cpu_time(end)
            print *, "... in ", end-start, " seconds."
        end if
!
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
        integer(kind=8) :: i_node, nb_nodes
!
        if (this%info >= 2) then
            print *, "Cleaning objects..."
        end if
!
        nb_nodes = size(this%nodes)
        do i_node = 1, nb_nodes
            if (allocated(this%nodes(i_node)%faces)) then
                deallocate (this%nodes(i_node)%faces)
            end if
            if (allocated(this%nodes(i_node)%edges)) then
                deallocate (this%nodes(i_node)%edges)
            end if
        end do
!
        deallocate (this%nodes)
        deallocate (this%edges)
        deallocate (this%faces)
        deallocate (this%volumes)
        deallocate (this%cells)
        deallocate (this%lastGhostsLayer)
        deallocate (this%faces_dege)
        deallocate (this%edges_dege)
!
    end subroutine
!
! ==================================================================================================
!
    subroutine numbering_edge(cell_type, nb_edge, edge_type, edge_loc)
!
        implicit none
!
        integer(kind=8), intent(in) :: cell_type
        integer(kind=8), intent(out) :: nb_edge, edge_type(12), edge_loc(3, 12)
! --------------------------------------------------------------------------------------------------
! Get edge connectivity of a cell
!
! --------------------------------------------------------------------------------------------------
        nb_edge = 0
        edge_type = 0
        edge_loc = 0
!
        if (cell_type == MT_QUAD4 .or. cell_type == MT_QUAD8 .or. cell_type == MT_QUAD9) then
            nb_edge = 4
            edge_loc(1:2, 1) = [1, 2]
            edge_loc(1:2, 2) = [2, 3]
            edge_loc(1:2, 3) = [3, 4]
            edge_loc(1:2, 4) = [4, 1]
            if (cell_type == MT_QUAD4) then
                edge_type(1:nb_edge) = MT_SEG2
            else
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3, 1:nb_edge) = [5, 6, 7, 8]
            end if
        elseif (cell_type == MT_TRIA3 .or. cell_type == MT_TRIA6 .or. cell_type == MT_TRIA7) then
            nb_edge = 3
            edge_loc(1:2, 1) = [1, 2]
            edge_loc(1:2, 2) = [2, 3]
            edge_loc(1:2, 3) = [3, 1]
            if (cell_type == MT_TRIA3) then
                edge_type(1:nb_edge) = MT_SEG2
            else
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3, 1:nb_edge) = [4, 5, 6]
            end if
        elseif (cell_type == MT_HEXA8 .or. cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
            nb_edge = 12
            edge_loc(1:2, 1) = [1, 2]
            edge_loc(1:2, 2) = [2, 3]
            edge_loc(1:2, 3) = [3, 4]
            edge_loc(1:2, 4) = [1, 4]
            edge_loc(1:2, 5) = [1, 5]
            edge_loc(1:2, 6) = [2, 6]
            edge_loc(1:2, 7) = [3, 7]
            edge_loc(1:2, 8) = [4, 8]
            edge_loc(1:2, 9) = [5, 6]
            edge_loc(1:2, 10) = [6, 7]
            edge_loc(1:2, 11) = [7, 8]
            edge_loc(1:2, 12) = [5, 8]
            if (cell_type == MT_HEXA8) then
                edge_type(1:nb_edge) = MT_SEG2
            elseif (cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3, 1:nb_edge) = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (cell_type == MT_TETRA4 .or. cell_type == MT_TETRA10 .or. &
                cell_type == MT_TETRA15) then
            nb_edge = 6
            edge_loc(1:2, 1) = [1, 2]
            edge_loc(1:2, 2) = [2, 3]
            edge_loc(1:2, 3) = [1, 3]
            edge_loc(1:2, 4) = [1, 4]
            edge_loc(1:2, 5) = [2, 4]
            edge_loc(1:2, 6) = [3, 4]
            if (cell_type == MT_TETRA4) then
                edge_type(1:nb_edge) = MT_SEG2
            elseif (cell_type == MT_TETRA10 .or. cell_type == MT_TETRA15) then
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3, 1:nb_edge) = [5, 6, 7, 8, 9, 10]
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (cell_type == MT_PENTA6 .or. cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 &
                .or. cell_type == MT_PENTA21) then
            nb_edge = 9
            edge_loc(1:2, 1) = [1, 2]
            edge_loc(1:2, 2) = [2, 3]
            edge_loc(1:2, 3) = [1, 3]
            edge_loc(1:2, 4) = [1, 4]
            edge_loc(1:2, 5) = [2, 5]
            edge_loc(1:2, 6) = [3, 6]
            edge_loc(1:2, 7) = [4, 5]
            edge_loc(1:2, 8) = [5, 6]
            edge_loc(1:2, 9) = [4, 6]
            if (cell_type == MT_PENTA6) then
                edge_type(1:nb_edge) = MT_SEG2
            elseif (cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 .or. &
                    cell_type == MT_PENTA21) then
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3, 1:nb_edge) = [7, 8, 9, 10, 11, 12, 13, 14, 15]
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (cell_type == MT_PYRAM5 .or. cell_type == MT_PYRAM13 .or. &
                cell_type == MT_PYRAM19) then
            nb_edge = 8
            edge_loc(1:2, 1) = [1, 2]
            edge_loc(1:2, 2) = [2, 3]
            edge_loc(1:2, 3) = [3, 4]
            edge_loc(1:2, 4) = [1, 4]
            edge_loc(1:2, 5) = [1, 5]
            edge_loc(1:2, 6) = [2, 5]
            edge_loc(1:2, 7) = [3, 5]
            edge_loc(1:2, 8) = [4, 5]
            if (cell_type == MT_PYRAM5) then
                edge_type(1:nb_edge) = MT_SEG2
            elseif (cell_type == MT_PYRAM13 .or. cell_type == MT_PYRAM19) then
                edge_type(1:nb_edge) = MT_SEG3
                edge_loc(3, 1:nb_edge) = [6, 7, 8, 9, 10, 11, 12, 13]
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
        integer(kind=8), intent(in) :: cell_type
        integer(kind=8), intent(out) :: nb_face, face_type(6), face_loc(9, 6)
! ---------------------------------------------------------------------------------
! Get face connectivity of a cell
!
! ---------------------------------------------------------------------------------
        nb_face = 0
        face_type = 0
        face_loc = 0
!
        if (cell_type == MT_HEXA8 .or. cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
            nb_face = 6
            face_loc(1:4, 1) = [1, 4, 3, 2]
            face_loc(1:4, 2) = [1, 2, 6, 5]
            face_loc(1:4, 3) = [2, 3, 7, 6]
            face_loc(1:4, 4) = [3, 4, 8, 7]
            face_loc(1:4, 5) = [1, 5, 8, 4]
            face_loc(1:4, 6) = [5, 6, 7, 8]
            if (cell_type == MT_HEXA8) then
                face_type(1:6) = MT_QUAD4
            elseif (cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
                face_type(1:6) = MT_QUAD8
                face_loc(5:8, 1) = [12, 11, 10, 9]
                face_loc(5:8, 2) = [9, 14, 17, 13]
                face_loc(5:8, 3) = [10, 15, 18, 14]
                face_loc(5:8, 4) = [11, 16, 19, 15]
                face_loc(5:8, 5) = [13, 20, 16, 12]
                face_loc(5:8, 6) = [17, 18, 19, 20]
                if (cell_type == MT_HEXA27) then
                    face_type(1:6) = MT_QUAD9
                    face_loc(9, 1:6) = [21, 22, 23, 24, 25, 26]
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (cell_type == MT_TETRA4 .or. cell_type == MT_TETRA10 &
                .or. cell_type == MT_TETRA15) then
            nb_face = 4
            face_loc(1:3, 1) = [1, 3, 2]
            face_loc(1:3, 2) = [1, 2, 4]
            face_loc(1:3, 3) = [1, 4, 3]
            face_loc(1:3, 4) = [2, 3, 4]
            if (cell_type == MT_TETRA4) then
                face_type(1:4) = MT_TRIA3
            elseif (cell_type == MT_TETRA10 .or. cell_type == MT_TETRA15) then
                face_type(1:4) = MT_TRIA6
                face_loc(4:6, 1) = [7, 6, 5]
                face_loc(4:6, 2) = [5, 9, 8]
                face_loc(4:6, 3) = [8, 10, 7]
                face_loc(4:6, 4) = [6, 10, 9]
                if (cell_type == MT_TETRA15) then
                    face_type(1:4) = MT_TRIA7
                    face_loc(7, 1:4) = [11, 12, 13, 14]
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (cell_type == MT_PENTA6 .or. cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 &
                .or. cell_type == MT_PENTA21) then
            nb_face = 5
            face_loc(1:4, 1) = [1, 2, 5, 4]
            face_loc(1:4, 2) = [2, 3, 6, 5]
            face_loc(1:4, 3) = [1, 4, 6, 3]
            face_loc(1:3, 4) = [1, 3, 2]
            face_loc(1:3, 5) = [4, 5, 6]
            if (cell_type == MT_PENTA6) then
                face_type(1:3) = MT_QUAD4
                face_type(4:5) = MT_TRIA3
            elseif (cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 .or. &
                    cell_type == MT_PENTA21) then
                face_type(1:3) = MT_QUAD8
                face_type(4:5) = MT_TRIA6
                face_loc(5:8, 1) = [7, 11, 13, 10]
                face_loc(5:8, 2) = [8, 12, 14, 11]
                face_loc(5:8, 3) = [10, 15, 12, 9]
                face_loc(4:6, 4) = [9, 8, 7]
                face_loc(4:6, 5) = [13, 14, 15]
                if (cell_type == MT_PENTA18 .or. cell_type == MT_PENTA21) then
                    face_type(1:3) = MT_QUAD9
                    face_loc(9, 1:3) = [16, 17, 18]
                    if (cell_type == MT_PENTA21) then
                        face_type(4:5) = MT_TRIA7
                        face_loc(7, 4:5) = [19, 20]
                    end if
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (cell_type == MT_PYRAM5 .or. cell_type == MT_PYRAM13 .or. &
                cell_type == MT_PYRAM19) then
            nb_face = 5
            face_loc(1:4, 1) = [1, 4, 3, 2]
            face_loc(1:3, 2) = [1, 2, 5]
            face_loc(1:3, 3) = [2, 3, 5]
            face_loc(1:3, 4) = [3, 4, 5]
            face_loc(1:3, 5) = [4, 1, 5]
            if (cell_type == MT_PYRAM5) then
                face_type(1) = MT_QUAD4
                face_type(2:5) = MT_TRIA3
            elseif (cell_type == MT_PYRAM13 .or. cell_type == MT_PYRAM19) then
                face_type(1) = MT_QUAD8
                face_type(2:5) = MT_TRIA6
                face_loc(5:8, 1) = [9, 8, 7, 6]
                face_loc(4:6, 2) = [6, 11, 10]
                face_loc(4:6, 3) = [7, 12, 11]
                face_loc(4:6, 4) = [8, 13, 12]
                face_loc(4:6, 5) = [9, 10, 13]
                if (cell_type == MT_PYRAM19) then
                    face_type(1) = MT_QUAD9
                    face_type(2:5) = MT_TRIA7
                    face_loc(9, 1) = 14
                    face_loc(7, 2:5) = [15, 16, 17, 18]
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
        integer(kind=8), intent(in) :: cell_type
        integer(kind=8), intent(out) :: nodes_loc(27)
! -----------------------------------------------------------------------------------
! Add special treatment
        integer(kind=8) :: i_node, nno
!
        nodes_loc = 0
        nno = this%converter%nno(cell_type)
!
        do i_node = 1, nno
            nodes_loc(i_node) = i_node
        end do
        if (cell_type .eq. MT_HEXA9) then
            nodes_loc(1:8) = [1, 2, 3, 4, 5, 6, 7, 8]
            nodes_loc(9) = 27
        end if
        if (cell_type .eq. MT_PENTA7) then
            nodes_loc(1:6) = [1, 2, 3, 4, 5, 6]
            nodes_loc(7) = 21
        end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine dividing_cell(cell_type, nb_sub, sub_type, sub_loc, conv_type)
!
        implicit none
!
        integer(kind=8), intent(in) :: cell_type
        integer(kind=8), intent(out) :: nb_sub, sub_type(10), sub_loc(10, 10), conv_type(10)
! --------------------------------------------------------------------------------------------------
! Get subdivision of a cell (node connectivity)
!
! --------------------------------------------------------------------------------------------------
        nb_sub = 0
        sub_type = 0
        sub_loc = 0
!
        if (cell_type == MT_SEG2 .or. cell_type == MT_SEG3) then
            nb_sub = 2
            sub_type(1:nb_sub) = MT_SEG2
            conv_type(1:nb_sub) = cell_type
            sub_loc(1:2, 1) = [1, 3]
            sub_loc(1:2, 2) = [2, 3]
        elseif (cell_type == MT_QUAD4 .or. cell_type == MT_QUAD8 .or. cell_type == MT_QUAD9) then
            nb_sub = 4
            sub_type(1:nb_sub) = MT_QUAD4
            conv_type(1:nb_sub) = cell_type
            sub_loc(1:4, 1) = [1, 5, 9, 8]
            sub_loc(1:4, 2) = [2, 6, 9, 5]
            sub_loc(1:4, 3) = [3, 7, 9, 6]
            sub_loc(1:4, 4) = [4, 8, 9, 7]
        elseif (cell_type == MT_TRIA3 .or. cell_type == MT_TRIA6 .or. cell_type == MT_TRIA7) then
            nb_sub = 4
            sub_type(1:nb_sub) = MT_TRIA3
            conv_type(1:nb_sub) = cell_type
            sub_loc(1:3, 1) = [1, 4, 6]
            sub_loc(1:3, 2) = [2, 5, 4]
            sub_loc(1:3, 3) = [3, 6, 5]
            sub_loc(1:3, 4) = [4, 5, 6]
        elseif (cell_type == MT_HEXA8 .or. cell_type == MT_HEXA20 .or. cell_type == MT_HEXA27) then
            nb_sub = 8
            sub_type(1:nb_sub) = MT_HEXA8
            conv_type(1:nb_sub) = cell_type
            sub_loc(1:8, 1) = [1, 9, 21, 12, 13, 22, 27, 25]
            sub_loc(1:8, 2) = [2, 10, 21, 9, 14, 23, 27, 22]
            sub_loc(1:8, 3) = [3, 11, 21, 10, 15, 24, 27, 23]
            sub_loc(1:8, 4) = [4, 12, 21, 11, 16, 25, 27, 24]
            sub_loc(1:8, 5) = [13, 22, 27, 25, 5, 17, 26, 20]
            sub_loc(1:8, 6) = [14, 23, 27, 22, 6, 18, 26, 17]
            sub_loc(1:8, 7) = [15, 24, 27, 23, 7, 19, 26, 18]
            sub_loc(1:8, 8) = [16, 25, 27, 24, 8, 20, 26, 19]
        elseif (cell_type == MT_TETRA4 .or. cell_type == MT_TETRA10 .or. &
                cell_type == MT_TETRA15) then
            nb_sub = 8
            conv_type(1:nb_sub) = cell_type
            sub_type(1:nb_sub) = MT_TETRA4
            sub_loc(1:4, 1) = [1, 5, 7, 8]
            sub_loc(1:4, 2) = [2, 6, 5, 9]
            sub_loc(1:4, 3) = [3, 10, 7, 6]
            sub_loc(1:4, 4) = [4, 8, 10, 9]
            sub_loc(1:4, 5) = [5, 6, 7, 9]
            sub_loc(1:4, 6) = [6, 9, 10, 7]
            sub_loc(1:4, 7) = [8, 9, 5, 7]
            sub_loc(1:4, 8) = [8, 10, 9, 7]
        elseif (cell_type == MT_PENTA6 .or. cell_type == MT_PENTA15 .or. cell_type == MT_PENTA18 &
                .or. cell_type == MT_PENTA21) then
            nb_sub = 8
            sub_type(1:nb_sub) = MT_PENTA6
            conv_type(1:nb_sub) = cell_type
            sub_loc(1:6, 1) = [1, 7, 9, 10, 16, 18]
            sub_loc(1:6, 2) = [7, 2, 8, 16, 11, 17]
            sub_loc(1:6, 3) = [9, 7, 8, 18, 16, 17]
            sub_loc(1:6, 4) = [9, 8, 3, 18, 17, 12]
            sub_loc(1:6, 5) = [10, 16, 18, 4, 13, 15]
            sub_loc(1:6, 6) = [16, 11, 17, 13, 5, 14]
            sub_loc(1:6, 7) = [18, 16, 17, 15, 13, 14]
            sub_loc(1:6, 8) = [18, 17, 12, 15, 14, 6]
        elseif (cell_type == MT_PYRAM5 .or. cell_type == MT_PYRAM13 .or. &
                cell_type == MT_PYRAM19) then
            nb_sub = 10
            sub_type(1:6) = MT_PYRAM5
            conv_type(1:6) = cell_type
            sub_loc(1:5, 1) = [1, 6, 14, 9, 10]
            sub_loc(1:5, 2) = [6, 2, 7, 14, 11]
            sub_loc(1:5, 3) = [3, 8, 14, 7, 12]
            sub_loc(1:5, 4) = [4, 9, 14, 8, 13]
            sub_loc(1:5, 5) = [10, 11, 12, 13, 5]
            sub_loc(1:5, 6) = [13, 12, 11, 10, 14]
            sub_type(7:10) = MT_TETRA4
            sub_loc(1:4, 7) = [6, 10, 11, 14]
            sub_loc(1:4, 8) = [7, 11, 12, 14]
            sub_loc(1:4, 9) = [8, 12, 13, 14]
            sub_loc(1:4, 10) = [9, 13, 10, 14]
            select case (cell_type)
            case (MT_PYRAM5)
                conv_type(7:10) = MT_TETRA4
            case (MT_PYRAM13)
                conv_type(7:10) = MT_TETRA10
            case (MT_PYRAM19)
                conv_type(7:10) = MT_TETRA15
            end select
        else
            ASSERT(ASTER_FALSE)
        end if

    end subroutine
!
! ==================================================================================================
!
    function owner_cell(this, nb_nodes, nodes)
!
        implicit none
!
        class(Mmesh), intent(in) :: this
        integer(kind=8), intent(in) :: nb_nodes
        integer(kind=8), intent(in) :: nodes(:)
        integer(kind=8) :: owner_cell
! -----------------------------------------------------------------------------------
! Add special treatment
        integer(kind=8) :: i_node
!
        owner_cell = this%nodes(nodes(1))%owner
!
        do i_node = 2, nb_nodes
            owner_cell = min(owner_cell, this%nodes(nodes(i_node))%owner)
        end do
!
    end function
!
! ==================================================================================================
!
    subroutine add_initial_cell(this, cell_id)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: cell_id
!
        integer(kind=8) :: cell_type, cell_dim, cell_nodes(27), nb_nodes, cell_index, i_node
!
        ASSERT(this%nb_total_cells < this%max_cells)
        this%nb_cells = this%nb_cells+1
        this%nb_total_cells = this%nb_total_cells+1
        cell_type = this%converter%map_type(this%v_typema(cell_id))
        cell_dim = this%converter%dim(cell_type)
!
        nb_nodes = this%converter%nno(cell_type)
        cell_nodes = 0
        do i_node = 1, nb_nodes
            cell_nodes(i_node) = this%connex_in(this%connex_map_in(cell_id)+i_node-1)
        end do
!
        if (this%debug) then
            print *, "Cell ", cell_id, ": ", cell_type, &
                this%converter%name(cell_type), cell_dim
        end if
!
        if (cell_dim == 3) then
            cell_index = this%add_volume(cell_type, cell_nodes)
            this%volumes(cell_index)%cell_id = cell_id
        elseif (cell_dim == 2) then
            cell_index = this%add_face(cell_type, cell_nodes)
            this%faces(cell_index)%cell_id = cell_id
        elseif (cell_dim == 1) then
            cell_index = this%add_edge(cell_type, cell_nodes)
            this%edges(cell_index)%cell_id = cell_id
        elseif (cell_dim == 0) then
            cell_index = this%add_point1(cell_type, cell_nodes)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        this%cells(cell_id)%type = cell_type
        this%cells(cell_id)%dim = cell_dim
        this%cells(cell_id)%id = cell_id
        this%cells(cell_id)%ss_id = cell_index
        this%cells(cell_id)%nodes(1:nb_nodes) = cell_nodes(1:nb_nodes)
!
    end subroutine
!
! ==================================================================================================
!
    function add_volume(this, type, nodes, parent, isub) result(volume_id)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: type, nodes(27)
        integer(kind=8), intent(in), optional :: parent, isub
        integer(kind=8) :: volume_id
! ----------------------------------------------------------------------
        integer(kind=8) :: nno, i_face
        integer(kind=8) :: nb_edges, edge_type(12), edge_loc(3, 12), edge_id, edge_nno
        integer(kind=8) :: nb_faces, face_type(6), face_loc(9, 6), i_node, i_edge
        integer(kind=8) :: face_nno, face_nodes(27), face_id, edge_nodes(27)
!
        ASSERT(this%converter%dim(type) == 3)
        nno = this%converter%nno(type)
!
        this%nb_volumes = this%nb_volumes+1
        if (this%nb_volumes > this%max_volumes) then
            call this%increase_memory("VOLUMES ", 2*this%max_volumes)
        end if
        ASSERT(this%nb_volumes <= this%max_volumes)
        volume_id = this%nb_volumes
        this%volumes(volume_id)%type = type
        this%volumes(volume_id)%nodes(1:nno) = nodes(1:nno)
        if (present(parent)) then
            this%volumes(volume_id)%parent = parent
        end if
        if (present(isub)) then
            this%volumes(volume_id)%isub = isub
        end if
! --- create edges
        call numbering_edge(type, nb_edges, edge_type, edge_loc)
        this%volumes(volume_id)%nb_edges = nb_edges
        do i_edge = 1, nb_edges
            edge_nno = this%converter%nno(edge_type(i_edge))
            edge_nodes = 0
            do i_node = 1, edge_nno
                edge_nodes(i_node) = nodes(edge_loc(i_node, i_edge))
            end do
            edge_id = this%add_edge(edge_type(i_edge), edge_nodes)
            this%volumes(volume_id)%edges(i_edge) = edge_id
            this%edges(edge_id)%nb_volumes = this%edges(edge_id)%nb_volumes+1
            this%edges(edge_id)%volumes(this%edges(edge_id)%nb_volumes) = volume_id
        end do
! --- create faces
        call numbering_face(type, nb_faces, face_type, face_loc)
        this%volumes(volume_id)%nb_faces = nb_faces
        do i_face = 1, nb_faces
            face_nno = this%converter%nno(face_type(i_face))
            face_nodes = 0
            do i_node = 1, face_nno
                face_nodes(i_node) = nodes(face_loc(i_node, i_face))
            end do
            face_id = this%add_face(face_type(i_face), face_nodes)
            this%volumes(volume_id)%faces(i_face) = face_id
            this%faces(face_id)%nb_volumes = this%faces(face_id)%nb_volumes+1
            this%faces(face_id)%volumes(this%faces(face_id)%nb_volumes) = volume_id
        end do
!
        call this%convert_volume(volume_id)
!
        if (this%debug) then
            print *, "Add volume: ", volume_id
            print *, "- Find: ", ASTER_FALSE
            print *, "- Type: ", this%converter%name(this%volumes(volume_id)%type), &
                "(", this%volumes(volume_id)%type, ")"
            print *, "- Owner: ", this%owner_cell(nno, this%volumes(volume_id)%nodes)
            print *, "- Nodes: ", this%volumes(volume_id)%nodes
            print *, "- Edges: ", this%volumes(volume_id)%edges
            print *, "- Faces: ", this%volumes(volume_id)%faces
            if (mult_elem(nno, this%volumes(volume_id)%nodes)) then
                ASSERT(ASTER_FALSE)
            end if
        end if
!
    end function
!
! ==================================================================================================
!
    function add_face(this, face_type, nodes, parent, isub) result(face_id)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: face_type, nodes(27)
        integer(kind=8), intent(in), optional :: parent, isub
        integer(kind=8) :: face_id
! ----------------------------------------------------------------------
        integer(kind=8) :: nno, nnos, nno_sort(9), i_edge, edge_id
        integer(kind=8) :: nb_edge, edge_type(12), edge_loc(3, 12), i_node
        integer(kind=8) :: edge_nno, edge_nodes(27), old_size
        integer(kind=8), allocatable :: new_faces(:)
        aster_logical :: find
!
        ASSERT(face_type > 0)
        ASSERT(this%converter%dim(face_type) == 2)
        nno = this%converter%nno(face_type)
        nnos = this%converter%nnos(face_type)
        nno_sort = 0
        nno_sort(1:nno) = nodes(1:nno)
!
        call sort_nodes_face(nno, nnos, nno_sort)
!
        call this%find_face(nno, nnos, nno_sort, find, face_id)
!
        if (.not. find) then
            this%nb_faces = this%nb_faces+1
            if (this%nb_faces > this%max_faces) then
                call this%increase_memory("FACES   ", 2*this%max_faces)
            end if
            ASSERT(this%nb_faces <= this%max_faces)
            face_id = this%nb_faces
            this%faces(face_id)%type = face_type
            this%faces(face_id)%nodes(1:nno) = nodes(1:nno)
            this%faces(face_id)%nnos = nnos
            this%faces(face_id)%nno_sort(1:nno) = nno_sort(1:nno)
            if (present(parent)) then
                this%faces(face_id)%parent = parent
            end if
            if (present(isub)) then
                this%faces(face_id)%isub = isub
            end if
! --- create edges
            call numbering_edge(face_type, nb_edge, edge_type, edge_loc)
            this%faces(face_id)%nb_edges = nb_edge
            do i_edge = 1, nb_edge
                edge_nno = this%converter%nno(edge_type(i_edge))
                do i_node = 1, edge_nno
                    edge_nodes(i_node) = nodes(edge_loc(i_node, i_edge))
                end do
                if (present(isub)) then
                    edge_id = this%add_edge(edge_type(i_edge), edge_nodes, &
                                            isub=i_edge, face_id=face_id)
                else
                    edge_id = this%add_edge(edge_type(i_edge), edge_nodes)
                end if
                this%faces(face_id)%edges(i_edge) = edge_id
                this%edges(edge_id)%nb_faces = this%edges(edge_id)%nb_faces+1
                ASSERT(this%edges(edge_id)%nb_faces <= 20)
                this%edges(edge_id)%faces(this%edges(edge_id)%nb_faces) = face_id
            end do
!
            if (this%nodes(nno_sort(1))%nb_faces >= this%nodes(nno_sort(1))%max_faces) then
                if (this%nodes(nno_sort(1))%max_faces > 0) then
                    old_size = this%nodes(nno_sort(1))%max_faces
                    allocate (new_faces(old_size))
                    new_faces(1:old_size) = this%nodes(nno_sort(1))%faces(1:old_size)
                    deallocate (this%nodes(nno_sort(1))%faces)
                    allocate (this%nodes(nno_sort(1))%faces(2*old_size))
                    this%nodes(nno_sort(1))%faces(1:old_size) = new_faces(1:old_size)
                    deallocate (new_faces)
                    this%nodes(nno_sort(1))%max_faces = 2*old_size
                else
                    this%nodes(nno_sort(1))%max_faces = 30
                    allocate (this%nodes(nno_sort(1))%faces(30))
                end if
            end if
            this%nodes(nno_sort(1))%nb_faces = this%nodes(nno_sort(1))%nb_faces+1
            this%nodes(nno_sort(1))%faces(this%nodes(nno_sort(1))%nb_faces) = face_id
!
            call this%convert_face(face_id)
        end if
!
        if (this%debug) then
            print *, "Add face: ", face_id
            print *, "- Find: ", find
            print *, "- Type: ", this%converter%name(this%faces(face_id)%type), &
                "(", this%faces(face_id)%type, ")"
            print *, "- Owner: ", this%owner_cell(nno, this%faces(face_id)%nodes)
            print *, "- Nodes: ", this%faces(face_id)%nodes
            print *, "- NNOS: ", this%faces(face_id)%nno_sort
            print *, "- Edges: ", this%faces(face_id)%edges
            if (mult_elem(nno, this%faces(face_id)%nodes)) then
                ASSERT(ASTER_FALSE)
            end if
        end if
!
    end function
!
! ==================================================================================================
!
    function add_edge(this, type, nodes, parent, isub, face_id) result(edge_id)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: type, nodes(27)
        integer(kind=8), intent(in), optional :: parent, isub, face_id
        integer(kind=8) :: edge_id
! ----------------------------------------------------------------------
        integer(kind=8) :: nno, nno_sort(3), old_size
        integer(kind=8), allocatable :: new_edges(:)
        aster_logical :: find
!
        ASSERT(this%converter%dim(type) == 1)
        nno = this%converter%nno(type)
        nno_sort(1) = minval(nodes(1:2))
        nno_sort(2) = maxval(nodes(1:2))
        if (nno > 2) then
            nno_sort(3:nno) = nodes(3:nno)
        end if
!
        call this%find_edge(nno, nno_sort, find, edge_id)
!
        if (.not. find) then
            this%nb_edges = this%nb_edges+1
            if (this%nb_edges > this%max_edges) then
                call this%increase_memory("EDGES   ", 2*this%max_edges)
            end if
            ASSERT(this%nb_edges <= this%max_edges)
            edge_id = this%nb_edges
            this%edges(edge_id)%type = type
            this%edges(edge_id)%nodes(1:nno) = nodes(1:nno)
            this%edges(edge_id)%nno_sort(1:nno) = nno_sort(1:nno)
            if (present(parent)) then
                this%edges(edge_id)%parent = parent
            end if
            if (present(isub)) then
                this%edges(edge_id)%isub = isub
            end if
!
            if (this%nodes(nno_sort(1))%nb_edges >= this%nodes(nno_sort(1))%max_edges) then
                if (this%nodes(nno_sort(1))%max_edges > 0) then
                    old_size = this%nodes(nno_sort(1))%max_edges
                    allocate (new_edges(old_size))
                    new_edges(1:old_size) = this%nodes(nno_sort(1))%edges(1:old_size)
                    deallocate (this%nodes(nno_sort(1))%edges)
                    allocate (this%nodes(nno_sort(1))%edges(2*old_size))
                    this%nodes(nno_sort(1))%edges(1:old_size) = new_edges(1:old_size)
                    deallocate (new_edges)
                    this%nodes(nno_sort(1))%max_edges = 2*old_size
                else
                    this%nodes(nno_sort(1))%max_edges = 30
                    allocate (this%nodes(nno_sort(1))%edges(30))
                end if
            end if
            this%nodes(nno_sort(1))%nb_edges = this%nodes(nno_sort(1))%nb_edges+1
            this%nodes(nno_sort(1))%edges(this%nodes(nno_sort(1))%nb_edges) = edge_id
!
            if (present(face_id)) then
                call this%convert_edge(edge_id, face_id)
            else
                call this%convert_edge(edge_id)
            end if
        end if
!
        if (this%debug) then
            print *, "Add edge: ", edge_id
            print *, "- Find: ", find
            print *, "- Type: ", this%converter%name(this%edges(edge_id)%type), &
                "(", this%edges(edge_id)%type, ")"
            print *, "- Owner: ", this%owner_cell(3, this%edges(edge_id)%nodes)
            print *, "- Nodes: ", this%edges(edge_id)%nodes(1:nno)
            print *, "- NNOS: ", this%edges(edge_id)%nno_sort(1:nno)
            if (mult_elem(nno, this%edges(edge_id)%nodes)) then
                ASSERT(ASTER_FALSE)
            end if

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
        integer(kind=8), intent(in) :: type, nodes(27)
        integer(kind=8) :: point_id
! ----------------------------------------------------------------------
!
        ASSERT(this%converter%dim(type) == 0)
        ASSERT(this%converter%nno(type) == 1)
        point_id = 1
!
        if (this%debug) then
            print *, "Add point1: "
            print *, "- Type: ", this%converter%name(type), "(", type, ")"
            print *, "- Nodes: ", nodes(1)
            print *, "- NNOS: ", nodes(1)
        end if
!
    end function
!
! ==================================================================================================
!
    function add_node(this, coor, owner) result(node_id)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        real(kind=8), intent(in) :: coor(3)
        integer(kind=8), intent(in) :: owner
        integer(kind=8) :: node_id
! ----------------------------------------------------------------------
!
        this%nb_nodes = this%nb_nodes+1
        this%nb_total_nodes = this%nb_total_nodes+1
        if (this%nb_nodes > this%max_nodes) then
            call this%increase_memory("NODES   ", 2*this%max_nodes)
            ASSERT(this%nb_nodes <= this%max_nodes)
        end if
        node_id = this%nb_nodes
        this%nodes(node_id)%id = node_id
        this%nodes(node_id)%keep = ASTER_TRUE
        this%nodes(node_id)%coor(1:3) = coor
        this%nodes(node_id)%owner = owner
!
    end function
!
! ==================================================================================================
!
    subroutine copy_mesh(this, mesh_out)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        character(len=8), intent(in) :: mesh_out
! ------------------------------------------------------------------
        character(len=24) :: cooval, coodsc, grpnoe
        character(len=24) :: gpptnn, grpmai, gpptnm, connex, titre, typmai, adapma
        character(len=4) :: dimesp
        integer(kind=8) :: i_node, nno, i_cell, ntgeo, nbnoma, node_id, nb_no_loc
        integer(kind=8) :: rank, nbproc, i_proc, deca, cell_id, pCellShift, hugeValue, iProc, iCount
        mpi_int :: mrank, msize
        real(kind=8):: start, end
        real(kind=8), pointer :: v_coor(:) => null()
        integer(kind=8), pointer :: v_int(:) => null()
        integer(kind=8), pointer :: v_connex(:) => null()
        integer(kind=8), pointer :: v_noex(:) => null()
        integer(kind=8), pointer :: v_maex(:) => null()
        integer(kind=8), pointer :: v_nuloc(:) => null()
        integer(kind=8), pointer :: v_nulogl(:) => null()
        integer(kind=8), pointer :: v_numagl(:) => null()
        integer(kind=8), pointer :: nbCellPerProc(:) => null()
        integer(kind=4), pointer :: v_lastlayer(:) => null()
!
        call jemarq()
!
        call this%check_mesh()
!
        call asmpi_info(rank=mrank, size=msize)
        rank = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
!
        if (this%info >= 2) then
            print *, "Copying mesh..."
            call cpu_time(start)
        end if
!
        call sdmail(mesh_out, cooval, coodsc, grpnoe, gpptnn, &
                    grpmai, gpptnm, connex, titre, typmai, &
                    adapma)
!
! --- Create nodes
!
! ------ Copy coordinates
        call wkvect(cooval, 'G V R', this%nb_nodes*3, vr=v_coor)
        if (this%isHPC) then
            call wkvect(mesh_out//'.NOEX', 'G V I', this%nb_nodes, vi=v_noex)
        end if
        call codent(this%dim_mesh, 'G', dimesp)
        call jeecra(cooval, 'DOCU', cval=dimesp)
        node_id = 0
        nb_no_loc = 0
        do i_node = 1, this%nb_total_nodes
            if (this%nodes(i_node)%keep) then
                node_id = node_id+1
                ASSERT(node_id == this%nodes(i_node)%id)
                v_coor(3*(node_id-1)+1:3*(node_id-1)+3) = this%nodes(i_node)%coor(1:3)
                if (this%isHPC) then
                    v_noex(node_id) = this%nodes(i_node)%owner
                    if (v_noex(node_id) == rank) then
                        nb_no_loc = nb_no_loc+1
                    end if
                end if
            end if
        end do
        ASSERT(node_id == this%nb_nodes)
! ------ Create Global numbering
        if (this%isHPC) then
            call wkvect('&&CREAMA.NOLOC', 'V V I', nbproc, vi=v_nuloc)
            v_nuloc(rank+1) = nb_no_loc
            call asmpi_comm_vect('MPI_SUM', 'I', nbval=nbproc, vi=v_nuloc)
            deca = 0
            do i_proc = 1, rank
                deca = deca+v_nuloc(i_proc)
            end do
!
            call wkvect(mesh_out//'.NUNOLG', 'G V I', this%nb_nodes, vi=v_nulogl)
            v_nulogl = -1
            do i_node = 1, this%nb_nodes
                if (v_noex(i_node) == rank) then
                    v_nulogl(i_node) = deca
                    deca = deca+1
                end if
            end do
            call jedetr('&&CREAMA.NOLOC')
        end if
! ------ Type of GEOM_R field
        call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), ntgeo)
        call wkvect(coodsc, 'G V I', 3, vi=v_int)
        call jeecra(coodsc, 'DOCU', 0, 'CHGO')
        v_int(1) = ntgeo
        v_int(2) = -3
        v_int(3) = 14
!
! --- Create cells
!
! ------ Count total number of nodes (repeated)
        nbnoma = 0
        do i_cell = 1, this%nb_total_cells
            if (this%cells(i_cell)%keep) then
                nbnoma = nbnoma+this%converter%nno(this%cells(i_cell)%type)
            end if
        end do
! ------ Create connectivity
        call wkvect(typmai, 'G V I', this%nb_cells, vi=v_int)
        call jecrec(connex, 'G V I', 'NU', 'CONTIG', 'VARIABLE', this%nb_cells)
        call jeecra(connex, 'NUTIOC', this%nb_cells)
        call jeecra(connex, 'LONT', nbnoma)
        if (this%isHPC) then
            call wkvect(mesh_out//'.MAEX', 'G V I', this%nb_cells, vi=v_maex)
            call wkvect(mesh_out//'.NUMALG', 'G V I', this%nb_cells, vi=v_numagl)
            AS_ALLOCATE(vi=nbCellPerProc, size=nbproc+1)
            nbCellPerProc(rank+2) = this%nb_total_cells
            call asmpi_comm_vect('MPI_SUM', 'I', nbval=nbproc+1, vi=nbCellPerProc)
            do iProc = 1, nbproc
                nbCellPerProc(iProc+1) = nbCellPerProc(iProc+1)+nbCellPerProc(iProc)
            end do
            pCellShift = nbCellPerProc(rank+1)
! --------- A huge value is used for global numbering of cells which are not
!           owned by current processor because to obtain the true value
!           it must be mandatory to communicate (#34152)
            hugeValue = huge(pCellShift)
            AS_DEALLOCATE(vi=nbCellPerProc)
        else
            pCellShift = 0
            hugeValue = 0
        end if
!
        iCount = pCellShift+1
        do i_cell = 1, this%nb_total_cells
            if (this%cells(i_cell)%keep) then
                cell_id = this%cells(i_cell)%id
                v_int(cell_id) = this%converter%cata_type(this%cells(i_cell)%type)
                nno = this%converter%nno(this%cells(i_cell)%type)
                call jeecra(jexnum(connex, cell_id), 'LONMAX', nno)
                call jeecra(jexnum(connex, cell_id), 'LONUTI', nno)
                call jeveuo(jexnum(connex, cell_id), 'E', vi=v_connex)
                do i_node = 1, nno
                    v_connex(i_node) = this%nodes(this%cells(i_cell)%nodes(i_node))%id
                end do
                if (this%isHPC) then
                    v_maex(cell_id) = this%owner_cell(nno, this%cells(i_cell)%nodes)
                    if (v_maex(cell_id) .eq. rank) then
                        v_numagl(cell_id) = iCount
                        iCount = iCount+1
                    else
                        v_numagl(cell_id) = hugeValue
                    end if
                end if
            end if
        end do
!
! --- Create groups
!
        call this%copy_group_no(grpnoe, gpptnn)
        call this%copy_group_ma(grpmai, gpptnm)
!
! --- Create .DIME
!
        call wkvect(mesh_out//'.DIME', 'G V I', 6, vi=v_int)
        v_int(1) = this%nb_nodes
        v_int(3) = this%nb_cells
        v_int(6) = this%dim_mesh
!
! --- Create joint for ParallelMesh
!
        call this%create_joints(mesh_out)
!
!
! --- Create .LASTGHOLAYER
!
        call wkvect(mesh_out//'.LASTGHOLAYER', 'G V S', max(1, this%lastLayerSize), vi4=v_lastlayer)
        call jeecra(mesh_out//'.LASTGHOLAYER', 'LONUTI', this%lastLayerSize)
        do i_node = 1, this%lastLayerSize
            v_lastlayer(i_node) = this%lastGhostsLayer(i_node)-1_4
        end do

        if (this%info >= 2) then
            call cpu_time(end)
            print *, "... in ", end-start, " seconds."
        end if
!
        call jedema()
!
    end subroutine
!
! ==================================================================================================
!
    subroutine convert_cells(this, nb_cells, list_cells)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: nb_cells, list_cells(nb_cells)
! ------------------------------------------------------------------
        integer(kind=8) :: i_cell, cell_id, cell_dim, object_id, nno, cell_type
        integer(kind=8) :: i_node, nodes_loc(27)
        real(kind=8):: start, end
!
        if (this%info >= 2) then
            print *, "Converting cells..."
            call cpu_time(start)
        end if
!
        do i_cell = 1, nb_cells
            cell_id = list_cells(i_cell)
            cell_dim = this%cells(cell_id)%dim
            cell_type = this%cells(cell_id)%type
            object_id = this%cells(cell_id)%ss_id
!
            if (this%debug) then
                print *, "Convert ", cell_id, ": ", this%cells(cell_id)%type, &
                    this%converter%name(this%cells(cell_id)%type), cell_dim
            end if
!
            if (this%converter%to_convert(cell_type)) then
                this%cells(cell_id)%type = this%converter%convert_to(cell_type)
                nno = this%converter%nno(this%cells(cell_id)%type)
                call this%numbering_nodes(this%cells(cell_id)%type, nodes_loc)
                this%cells(cell_id)%nodes = 0
                if (cell_dim == 3) then
                    do i_node = 1, nno
                        this%cells(cell_id)%nodes(i_node) = &
                            this%volumes(object_id)%nodes(nodes_loc(i_node))
                    end do
                elseif (cell_dim == 2) then
                    do i_node = 1, nno
                        this%cells(cell_id)%nodes(i_node) = &
                            this%faces(object_id)%nodes(nodes_loc(i_node))
                    end do
                elseif (cell_dim == 1) then
                    do i_node = 1, nno
                        this%cells(cell_id)%nodes(i_node) = &
                            this%edges(object_id)%nodes(nodes_loc(i_node))
                    end do
                elseif (cell_dim == 0) then
                    ASSERT(ASTER_FALSE)
                else
                    ASSERT(ASTER_FALSE)
                end if
            end if
        end do
! --- Keep only necessary nodes
        call this%update()
!
        if (this%info >= 2) then
            call cpu_time(end)
            print *, "... in ", end-start, " seconds."
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
        integer(kind=8), intent(in) :: volu_id
! ------------------------------------------------------------------
        integer(kind=8) :: volu_type, volu_type_end
        integer(kind=8) :: nno, nno_end, node_id, i_face, i_node, face_type, face_nno
        integer(kind=8) :: nb_edges, edge_type(12), edge_loc(3, 12), face_id
        integer(kind=8) :: nb_faces, faces_type(6), face_loc(9, 6), i_edge, owner
!
        volu_type = this%volumes(volu_id)%type
        nno = this%converter%nno(volu_type)
        volu_type_end = this%converter%convert_max(volu_type)
        nno_end = this%converter%nno(volu_type_end)
!
        if (this%debug) then
            print *, "Volume: ", volu_id, volu_type, volu_type_end
        end if
!
        if (nno_end > nno) then
! --- Add nodes on edges
            if (volu_type == MT_HEXA8 .or. volu_type == MT_TETRA4 .or. &
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
            owner = this%owner_cell(nno_end-1, this%volumes(volu_id)%nodes)
            node_id = this%add_node(this%barycenter(volu_id, 3), owner)
            this%volumes(volu_id)%nodes(nno_end) = node_id
        else
            ASSERT(volu_type == volu_type_end)
        end if
!
        this%volumes(volu_id)%type = volu_type_end
!
        if (this%debug) then
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
        integer(kind=8), intent(in) :: face_id
! ------------------------------------------------------------------
        integer(kind=8) :: face_type, face_type_end, owner
        integer(kind=8) :: nno, nno_end, node_id, i_edge, i_node, nnos
!
        face_type = this%faces(face_id)%type
        nno = this%converter%nno(face_type)
        nnos = this%converter%nnos(face_type)
        face_type_end = this%converter%convert_max(face_type)
        nno_end = this%converter%nno(face_type_end)
!
        if (this%debug) then
            print *, "Face: ", face_id, face_type, face_type_end, this%faces(face_id)%nodes
        end if
!
        if (nno_end > nno) then
! --- Add nodes at the middle of edges
            if (face_type == MT_TRIA3 .or. face_type == MT_QUAD4) then
                do i_edge = 1, this%faces(face_id)%nb_edges
                    call this%convert_edge(this%faces(face_id)%edges(i_edge), face_id)
                    this%faces(face_id)%nodes(nno+i_edge) = &
                        this%edges(this%faces(face_id)%edges(i_edge))%nodes(3)
                end do
                this%faces(face_id)%nno_sort(nnos+1:2*nnos) = &
                    this%faces(face_id)%nodes(nnos+1:2*nnos)
                call qsort(this%faces(face_id)%nno_sort(nnos+1:2*nnos))
            end if
! --- Add node at the barycenter
            owner = this%owner_cell(nno_end-1, this%faces(face_id)%nodes)
            node_id = this%add_node(this%barycenter(face_id, 2), owner)
            this%faces(face_id)%nodes(nno_end) = node_id
        else
            ASSERT(face_type == face_type_end)
        end if
!
        this%faces(face_id)%type = face_type_end
!
        if (this%debug) then
            do i_node = 1, nno_end
                ASSERT(this%faces(face_id)%nodes(i_node) > 0)
            end do
        end if
!
    end subroutine
!
! ==================================================================================================
!
    subroutine convert_edge(this, edge_id, face_id)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: edge_id
        integer(kind=8), intent(in), optional :: face_id
! ------------------------------------------------------------------
        integer(kind=8) :: edge_type, edge_type_end, owner
        integer(kind=8) :: nno, nno_end, node_id, i_node
!
        edge_type = this%edges(edge_id)%type
        nno = this%converter%nno(edge_type)
        edge_type_end = this%converter%convert_max(edge_type)
        nno_end = this%converter%nno(edge_type_end)
!
        if (this%debug) then
            print *, "Edge: ", edge_id, edge_type, edge_type_end, this%edges(edge_id)%nodes
        end if
!
        if (nno_end > nno) then
            if (edge_type_end == MT_SEG3) then
                owner = this%owner_cell(2, this%edges(edge_id)%nodes)
                if (present(face_id)) then
                    node_id = this%add_node(this%barycenter_fe(edge_id, face_id), owner)
                else
                    node_id = this%add_node(this%barycenter(edge_id, 1), owner)
                end if
                this%edges(edge_id)%nodes(3) = node_id
                this%edges(edge_id)%nno_sort(3) = node_id
            else
                ASSERT(ASTER_FALSE)
                node_id = this%add_node([0.d0, 0.d0, 0.d0], -1)
                this%edges(edge_id)%nodes(3) = node_id
                this%edges(edge_id)%nno_sort(3) = node_id
                node_id = this%add_node([0.d0, 0.d0, 0.d0], -1)
                this%edges(edge_id)%nodes(4) = node_id
                this%edges(edge_id)%nno_sort(4) = node_id
                call qsort(this%edges(edge_id)%nno_sort(3:4))
            end if
        else
            ASSERT(edge_type == edge_type_end)
        end if
!
        this%edges(edge_id)%type = edge_type_end
!
        if (this%debug) then
            do i_node = 1, nno_end
                ASSERT(this%edges(edge_id)%nodes(i_node) > 0)
            end do
        end if
!
    end subroutine
!
! ==================================================================================================
!
    function barycenter(this, cell_id, cell_dim) result(coor)
!
        implicit none
!
        class(Mmesh), intent(in) :: this
        integer(kind=8), intent(in) :: cell_id, cell_dim
        real(kind=8) :: coor(3)
! ---------------------------------------------------------------------------------
!
        if (cell_dim == 3) then
            coor = this%barycenter_volume(cell_id)
        elseif (cell_dim == 2) then
            coor = this%barycenter_face(cell_id)
        elseif (cell_dim == 1) then
            coor = this%barycenter_edge(cell_id)
        else
            ASSERT(ASTER_FALSE)
        end if
    end function
!
! ==================================================================================================
!
    function barycenter_edge(this, edge_id) result(coor)
!
        implicit none
!
        class(Mmesh), intent(in) :: this
        integer(kind=8), intent(in) :: edge_id
        real(kind=8) :: coor(3)
! ---------------------------------------------------------------------------------
        integer(kind=8) :: i_node, node
        real(kind=8) :: basis(27)
!
        coor = 0.d0
!
        if (this%edges(edge_id)%isub == 0) then
            ! coor_ref = [0.d0, 0.d0, 0.d0]
            ! call elrfvf("SE2", coor_ref, basis)
            basis(1:2) = 0.5d0
            do i_node = 1, 2
                node = this%edges(edge_id)%nodes(i_node)
                coor(1:3) = coor(1:3)+this%nodes(node)%coor(1:3)*basis(i_node)
            end do
        else
            if (this%edges(edge_id)%isub == 1) then
                ! coor_ref = [-0.5d0, 0.d0, 0.d0]
                ! call elrfvf("SE3", coor_ref, basis)
                basis(1:3) = [0.375d0, -0.125d0, 0.75d0]
            else
                ! coor_ref = [0.5d0, 0.d0, 0.d0]
                ! call elrfvf("SE3", coor_ref, basis)
                basis(1:3) = [-0.125, 0.375d0, 0.75d0]
            end if
            ASSERT(this%edges(edge_id)%parent > 0)
            do i_node = 1, 3
                node = this%edges(this%edges(edge_id)%parent)%nodes(i_node)
                coor(1:3) = coor(1:3)+this%nodes(node)%coor(1:3)*basis(i_node)
            end do
        end if
!
    end function
!
! ==================================================================================================
!
    function barycenter_fe(this, edge_id, face_id) result(coor)
!
! --- Compute barycenter of a edge using geometry of the face
!     It allows to preserve curvature
!
        implicit none
!
        class(Mmesh), intent(in) :: this
        integer(kind=8), intent(in) :: edge_id, face_id
        real(kind=8) :: coor(3)
! ---------------------------------------------------------------------------------
        integer(kind=8) :: i_node, node, type, nbnode, parent
        real(kind=8) :: basis(27), coor_ref(3)
!
        coor = 0.d0
!
        type = this%faces(face_id)%type
        if (type == MT_TRIA3 .or. type == MT_TRIA6 .or. type == MT_TRIA7) then
            if (this%faces(face_id)%isub == 1) then
                if (this%edges(edge_id)%isub == 1) then
                    coor_ref = [0.25d0, 0.d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 2) then
                    coor_ref = [0.25d0, 0.25d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 3) then
                    coor_ref = [0.d0, 0.25d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (this%faces(face_id)%isub == 2) then
                if (this%edges(edge_id)%isub == 1) then
                    coor_ref = [0.75d0, 0.25d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 2) then
                    coor_ref = [0.5d0, 0.25d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 3) then
                    coor_ref = [0.75d0, 0.d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (this%faces(face_id)%isub == 3) then
                if (this%edges(edge_id)%isub == 1) then
                    coor_ref = [0.d0, 0.75d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 2) then
                    coor_ref = [0.25d0, 0.5d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 3) then
                    coor_ref = [0.25d0, 0.75d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (this%faces(face_id)%isub == 4) then
                if (this%edges(edge_id)%isub == 1) then
                    coor_ref = [0.5d0, 0.25d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 2) then
                    coor_ref = [0.25d0, 0.5d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 3) then
                    coor_ref = [0.25d0, 0.25d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        elseif (type == MT_QUAD4 .or. type == MT_QUAD8 .or. type == MT_QUAD9) then
            if (this%faces(face_id)%isub == 1) then
                if (this%edges(edge_id)%isub == 1) then
                    coor_ref = [-0.5d0, -1.d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 2) then
                    coor_ref = [0.0d0, -0.5d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 3) then
                    coor_ref = [-0.5d0, 0.d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 4) then
                    coor_ref = [-1.d0, -0.5d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (this%faces(face_id)%isub == 2) then
                if (this%edges(edge_id)%isub == 1) then
                    coor_ref = [1.0d0, -0.5d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 2) then
                    coor_ref = [0.5d0, 0.0d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 3) then
                    coor_ref = [0.0d0, -0.5d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 4) then
                    coor_ref = [0.5d0, -1.0d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (this%faces(face_id)%isub == 3) then
                if (this%edges(edge_id)%isub == 1) then
                    coor_ref = [0.5d0, 1.d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 2) then
                    coor_ref = [0.0d0, 0.5d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 3) then
                    coor_ref = [0.5d0, 0.d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 4) then
                    coor_ref = [1.d0, 0.5d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (this%faces(face_id)%isub == 4) then
                if (this%edges(edge_id)%isub == 1) then
                    coor_ref = [-1.0d0, 0.5d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 2) then
                    coor_ref = [-0.5d0, 0.d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 3) then
                    coor_ref = [0.d0, 0.5d0, 0.d0]
                elseif (this%edges(edge_id)%isub == 4) then
                    coor_ref = [-0.5d0, 1.0d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

        parent = this%faces(face_id)%parent
        ASSERT(parent > 0)
        call elrfvf(this%converter%short_name(this%faces(parent)%type), coor_ref, basis, nbnode)
        do i_node = 1, nbnode
            node = this%faces(parent)%nodes(i_node)
            coor(1:3) = coor(1:3)+this%nodes(node)%coor(1:3)*basis(i_node)
        end do
!
    end function
!
! ==================================================================================================
!
    function barycenter_face(this, face_id) result(coor)
!
        implicit none
!
        class(Mmesh), intent(in) :: this
        integer(kind=8), intent(in) :: face_id
        real(kind=8) :: coor(3)
! ---------------------------------------------------------------------------------
        integer(kind=8) :: i_node, node, nbnode, parent, type
        real(kind=8) :: basis(27), coor_ref(3)
        character(len=8) :: stype
!
        coor = 0.d0
!
        type = this%faces(face_id)%type
        stype = this%converter%short_name(type)
        if (this%faces(face_id)%isub == 0) then
            if (type == MT_TRIA3 .or. type == MT_TRIA6 .or. type == MT_TRIA7) then
                coor_ref = [1.d0/3.d0, 1.d0/3.d0, 0.d0]
            elseif (type == MT_QUAD4 .or. type == MT_QUAD8 .or. type == MT_QUAD9) then
                coor_ref = [0.d0, 0.d0, 0.d0]
            else
                ASSERT(ASTER_FALSE)
            end if
            call elrfvf(stype, coor_ref, basis, nbnode)
            do i_node = 1, nbnode
                node = this%faces(face_id)%nodes(i_node)
                coor(1:3) = coor(1:3)+this%nodes(node)%coor(1:3)*basis(i_node)
            end do
        else
            if (type == MT_TRIA3 .or. type == MT_TRIA6 .or. type == MT_TRIA7) then
                if (this%faces(face_id)%isub == 1) then
                    coor_ref = [0.5d0/3.d0, 0.5d0/3.d0, 0.d0]
                elseif (this%faces(face_id)%isub == 2) then
                    coor_ref = [2.d0/3.d0, 0.5d0/3.d0, 0.d0]
                elseif (this%faces(face_id)%isub == 3) then
                    coor_ref = [0.5d0/3.d0, 2.d0/3.d0, 0.d0]
                elseif (this%faces(face_id)%isub == 4) then
                    coor_ref = [1.d0/3.d0, 1.d0/3.d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (type == MT_QUAD4 .or. type == MT_QUAD8 .or. type == MT_QUAD9) then
                if (this%faces(face_id)%isub == 1) then
                    coor_ref = [-0.5d0, -0.5d0, 0.d0]
                elseif (this%faces(face_id)%isub == 2) then
                    coor_ref = [0.5d0, -0.5d0, 0.d0]
                elseif (this%faces(face_id)%isub == 3) then
                    coor_ref = [0.5d0, 0.5d0, 0.d0]
                elseif (this%faces(face_id)%isub == 4) then
                    coor_ref = [-0.5d0, 0.5d0, 0.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            else
                ASSERT(ASTER_FALSE)
            end if

            parent = this%faces(face_id)%parent
            ASSERT(parent > 0)
            call elrfvf(this%converter%short_name(this%faces(parent)%type), coor_ref, basis, &
                        nbnode)
            do i_node = 1, nbnode
                node = this%faces(parent)%nodes(i_node)
                coor(1:3) = coor(1:3)+this%nodes(node)%coor(1:3)*basis(i_node)
            end do
        end if
!
    end function
!
! ==================================================================================================
!
    function barycenter_volume(this, volu_id) result(coor)
!
        implicit none
!
        class(Mmesh), intent(in) :: this
        integer(kind=8), intent(in) :: volu_id
        real(kind=8) :: coor(3)
! ---------------------------------------------------------------------------------
        integer(kind=8) :: i_node, node, nbnode, parent, type
        real(kind=8) :: basis(27), coor_ref(3)
        character(len=8) :: stype
!
        coor = 0.d0
!
        type = this%volumes(volu_id)%type
        stype = this%converter%short_name(type)
        if (this%volumes(volu_id)%isub == 0) then
            if (type == MT_TETRA4 .or. type == MT_TETRA10 .or. type == MT_TETRA15) then
                coor_ref = [0.25d0, 0.25d0, 0.25d0]
            elseif (type == MT_HEXA8 .or. type == MT_HEXA20 .or. type == MT_HEXA27) then
                coor_ref = [0.d0, 0.d0, 0.d0]
            elseif (type == MT_PENTA6 .or. type == MT_PENTA15 .or. type == MT_PENTA18 &
                    .or. type == MT_PENTA21) then
                coor_ref = [0.d0, 1.d0/3.d0, 1.d0/3.d0]
            elseif (type == MT_PYRAM5 .or. type == MT_PYRAM13 .or. type == MT_PYRAM19) then
                coor_ref = [0.d0, 0.0, 0.2d0]
            else
                ASSERT(ASTER_FALSE)
            end if
            call elrfvf(stype, coor_ref, basis, nbnode)
            do i_node = 1, nbnode
                node = this%volumes(volu_id)%nodes(i_node)
                coor(1:3) = coor(1:3)+this%nodes(node)%coor(1:3)*basis(i_node)
            end do
        else
            if (type == MT_TETRA4 .or. type == MT_TETRA10 .or. type == MT_TETRA15) then
                if (this%volumes(volu_id)%isub == 1) then
                    coor_ref = [0.5d0/4.d0, 2.5d0/4.d0, 0.5d0/4.d0]
                elseif (this%volumes(volu_id)%isub == 2) then
                    coor_ref = [0.5d0/4.d0, 0.5d0/4.d0, 2.5d0/4.d0]
                elseif (this%volumes(volu_id)%isub == 3) then
                    coor_ref = [0.5d0/4.d0, 0.5d0/4.d0, 0.5d0/4.d0]
                elseif (this%volumes(volu_id)%isub == 4) then
                    coor_ref = [2.5d0/4.d0, 0.5d0/4.d0, 0.5d0/4.d0]
                elseif (this%volumes(volu_id)%isub == 5) then
                    coor_ref = [0.5d0/4.d0, 1.d0/4.d0, 1.5d0/4.d0]
                elseif (this%volumes(volu_id)%isub == 6) then
                    coor_ref = [1.d0/4.d0, 0.5d0/4.d0, 1.d0/4.d0]
                elseif (this%volumes(volu_id)%isub == 7) then
                    coor_ref = [1.d0/4.d0, 1.5d0/4.d0, 1.d0/4.d0]
                elseif (this%volumes(volu_id)%isub == 8) then
                    coor_ref = [1.5d0/4.d0, 1.d0/4.d0, 0.5d0/4.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (type == MT_HEXA8 .or. type == MT_HEXA20 .or. type == MT_HEXA27) then
                if (this%volumes(volu_id)%isub == 1) then
                    coor_ref = [-0.5d0, -0.5d0, -0.5d0]
                elseif (this%volumes(volu_id)%isub == 2) then
                    coor_ref = [0.5d0, -0.5d0, -0.5d0]
                elseif (this%volumes(volu_id)%isub == 3) then
                    coor_ref = [0.5d0, 0.5d0, -0.5d0]
                elseif (this%volumes(volu_id)%isub == 4) then
                    coor_ref = [-0.5d0, 0.5d0, -0.5d0]
                elseif (this%volumes(volu_id)%isub == 5) then
                    coor_ref = [-0.5d0, -0.5d0, 0.5d0]
                elseif (this%volumes(volu_id)%isub == 6) then
                    coor_ref = [0.5d0, -0.5d0, 0.5d0]
                elseif (this%volumes(volu_id)%isub == 7) then
                    coor_ref = [0.5d0, 0.5d0, 0.5d0]
                elseif (this%volumes(volu_id)%isub == 8) then
                    coor_ref = [-0.5d0, 0.5d0, 0.5d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (type == MT_PENTA6 .or. type == MT_PENTA15 .or. type == MT_PENTA18 &
                    .or. type == MT_PENTA21) then
                if (this%volumes(volu_id)%isub == 1) then
                    coor_ref = [-0.5d0, 2.d0/3.d0, 0.5d0/3.d0]
                elseif (this%volumes(volu_id)%isub == 2) then
                    coor_ref = [-0.5d0, 0.5d0/3.d0, 2.d0/3.d0]
                elseif (this%volumes(volu_id)%isub == 3) then
                    coor_ref = [-0.5d0, 1.d0/3.d0, 1.d0/3.d0]
                elseif (this%volumes(volu_id)%isub == 4) then
                    coor_ref = [-0.5d0, 0.5d0/3.d0, 0.5d0/3.d0]
                elseif (this%volumes(volu_id)%isub == 5) then
                    coor_ref = [0.5d0, 2.d0/3.d0, 0.5d0/3.d0]
                elseif (this%volumes(volu_id)%isub == 6) then
                    coor_ref = [0.5d0, 0.5d0/3.d0, 2.d0/3.d0]
                elseif (this%volumes(volu_id)%isub == 7) then
                    coor_ref = [0.5d0, 1.d0/3.d0, 1.d0/3.d0]
                elseif (this%volumes(volu_id)%isub == 8) then
                    coor_ref = [0.5d0, 0.5d0/3.d0, 0.5d0/3.d0]
                else
                    ASSERT(ASTER_FALSE)
                end if
            elseif (type == MT_PYRAM5 .or. type == MT_PYRAM13 .or. type == MT_PYRAM19) then
                ASSERT(ASTER_FALSE)
                ! if (this%volumes(volu_id)%isub == 1) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 2) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 3) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 4) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 5) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 6) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 7) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 8) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 9) then
                !     coor_ref = []
                ! elseif (this%volumes(volu_id)%isub == 10) then
                !     coor_ref = []
                ! else
                !     ASSERT(ASTER_FALSE)
                ! end if
            else
                ASSERT(ASTER_FALSE)
            end if

            parent = this%volumes(volu_id)%parent
            ASSERT(parent > 0)
            call elrfvf(this%converter%short_name(this%volumes(parent)%type), &
                        coor_ref, basis, nbnode)
            do i_node = 1, nbnode
                node = this%volumes(parent)%nodes(i_node)
                coor(1:3) = coor(1:3)+this%nodes(node)%coor(1:3)*basis(i_node)
            end do
        end if
!
    end function
!
! ==================================================================================================
!
    subroutine update(this)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
! -----------------------------------------------------------------------
        integer(kind=8) :: i_node, i_cell, nno, node_id, rank, i
        integer(kind=8) :: i_layer, nb_ghost_cell, i_ghost_cell
        integer(kind=8), pointer :: ghost_cells(:) => null()
        integer(kind=8), pointer :: old2new(:) => null()
        mpi_int :: mrank
        aster_logical :: keep, already_exists
!
        call asmpi_info(rank=mrank)
        rank = to_aster_int(mrank)
!
! --- Keep initial orphelan nodes
        do i_node = 1, this%nb_total_nodes
            this%nodes(i_node)%keep = ASTER_FALSE
            if (this%nodes(i_node)%orphelan) then
                this%nodes(i_node)%keep = ASTER_TRUE
            end if
        end do
!
! --- Remove ghost cells - added later
        nb_ghost_cell = 0
        if (this%isHPC) then
            AS_ALLOCATE(vi=ghost_cells, size=this%nb_total_cells)
            do i_cell = 1, this%nb_total_cells
                if (this%cells(i_cell)%keep) then
                    nno = this%converter%nno(this%cells(i_cell)%type)
!
                    keep = ASTER_TRUE
                    do i_node = 1, nno
                        node_id = this%cells(i_cell)%nodes(i_node)
                        if (this%nodes(node_id)%owner /= rank) then
                            keep = ASTER_FALSE
                            nb_ghost_cell = nb_ghost_cell+1
                            ghost_cells(nb_ghost_cell) = i_cell
                            exit
                        end if
                    end do
                    this%cells(i_cell)%keep = keep
                end if
            end do
        end if
!
! --- Keep only nodes of cells
        do i_cell = 1, this%nb_total_cells
            if (this%cells(i_cell)%keep) then
                nno = this%converter%nno(this%cells(i_cell)%type)
!
                if (this%debug) then
                    print *, "Cell: ", i_cell, this%cells(i_cell)%type, nno, &
                        this%cells(i_cell)%nodes(1:nno)
                end if
!
                do i_node = 1, nno
                    node_id = this%cells(i_cell)%nodes(i_node)
                    this%nodes(node_id)%keep = ASTER_TRUE
                end do
            end if
        end do
!
        if (this%isHPC) then
! --- Add ghost layers
            do i_layer = 1, this%nb_layer
                do i_ghost_cell = 1, nb_ghost_cell
                    i_cell = ghost_cells(i_ghost_cell)
                    if (i_cell > 0) then
                        if (.not. this%cells(i_cell)%keep) then
                            nno = this%converter%nno(this%cells(i_cell)%type)
!
                            do i_node = 1, nno
                                node_id = this%cells(i_cell)%nodes(i_node)

                                if (this%nodes(node_id)%keep .or. &
                                    this%nodes(node_id)%owner == rank) then
                                    this%cells(i_cell)%keep = ASTER_TRUE
                                    exit
                                end if
                            end do
                        end if
                    end if
                end do
!
! --- Keep only nodes of cells
                this%lastLayerSize = 0
                do i_ghost_cell = 1, nb_ghost_cell
                    i_cell = ghost_cells(i_ghost_cell)
                    if (i_cell > 0) then
                        if (this%cells(i_cell)%keep) then
                            nno = this%converter%nno(this%cells(i_cell)%type)
!
                            if (this%debug) then
                                print *, "Cell: ", i_cell, this%cells(i_cell)%type, nno, &
                                    this%cells(i_cell)%nodes(1:nno)
                            end if
!
                            do i_node = 1, nno
                                node_id = this%cells(i_cell)%nodes(i_node)
                                if (i_layer == this%nb_layer .and. &
                                    .not. this%nodes(node_id)%keep) then
                                    ! Check if the value already exists
                                    already_exists = .false.
                                    do i = 1, this%lastLayerSize
                                        if (this%lastGhostsLayer(i) == &
                                            int(this%nodes(node_id)%id, 4)) then
                                            already_exists = .true.
                                            cycle
                                        end if
                                    end do

                                    ! Insert only if it doesn't exist
                                    if (.not. already_exists) then
                                        this%lastLayerSize = this%lastLayerSize+1
                                        this%lastGhostsLayer(this%lastLayerSize) = &
                                            int(this%nodes(node_id)%id, 4)
                                    end if
                                end if
                                this%nodes(node_id)%keep = ASTER_TRUE
                            end do
                            ghost_cells(i_ghost_cell) = -1
                        end if
                    end if
                end do
            end do
!
            AS_DEALLOCATE(vi=ghost_cells)
        end if
!
! --- Renumbering nodes
        this%nb_nodes = 0
        AS_ALLOCATE(vi=old2new, size=this%nb_total_nodes)
        do i_node = 1, this%nb_total_nodes
            if (this%nodes(i_node)%keep) then
                this%nb_nodes = this%nb_nodes+1
                old2new(this%nodes(i_node)%id) = this%nb_nodes
                this%nodes(i_node)%id = this%nb_nodes
            end if
        end do
!
! --- Renumbering last layer
        if (this%isHPC) then
            do i_node = 1, this%lastLayerSize
                this%lastGhostsLayer(i_node) = int(old2new(this%lastGhostsLayer(i_node)), 4)
            end do
        end if
        AS_DEALLOCATE(vi=old2new)
!
! --- Renumbering cells
        this%nb_cells = 0
        do i_cell = 1, this%nb_total_cells
            if (this%cells(i_cell)%keep) then
                this%nb_cells = this%nb_cells+1
                this%cells(i_cell)%id = this%nb_cells
            end if
        end do
!
        if (this%debug) then
            print *, "Update nodes: ", this%nb_nodes, this%nb_total_nodes
            do i_node = 1, this%nb_total_nodes
                if (this%nodes(i_node)%keep) then
                    print *, "Node: ", i_node, this%nodes(i_node)%id, this%nodes(i_node)%owner, &
                        this%nodes(i_node)%coor
                end if
            end do
            print *, "Update cells: ", this%nb_cells, this%nb_total_cells
            do i_cell = 1, this%nb_total_cells
                if (this%cells(i_cell)%keep) then
                    print *, "Cell: ", i_cell, this%cells(i_cell)%id, this%cells(i_cell)%type
                end if
            end do
        end if
        ASSERT(this%nb_nodes <= this%nb_total_nodes)
        ASSERT(this%nb_nodes > 0)
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
        integer(kind=8) :: i_cell, i_node, i_edge, i_face, i_volume, nno, nno1, nno2
        integer(kind=8) :: nb_edges, edges_type(12), edges_loc(3, 12), edge_id
        integer(kind=8) :: face_type, volu_type, face_id
        integer(kind=8) :: nb_faces, faces_type(6), faces_loc(9, 6)
!
! --- Check Nodes
        do i_cell = 1, this%nb_total_cells
            if (this%cells(i_cell)%keep) then
                nno = this%converter%nno(this%cells(i_cell)%type)
                do i_node = 1, nno
                    ASSERT(this%nodes(this%cells(i_cell)%nodes(i_node))%keep)
                    ASSERT(this%nodes(this%cells(i_cell)%nodes(i_node))%owner >= 0)
                end do
            end if
        end do
!
! --- Check Edges
        do i_edge = 1, this%nb_edges
            nno = this%converter%nno(this%edges(i_edge)%type)
            do i_node = 1, nno
                ASSERT(this%edges(i_edge)%nodes(i_node) > 0)
            end do
!
            ASSERT(this%edges(i_edge)%nb_volumes <= 25)
            do i_volume = 1, this%edges(i_edge)%nb_volumes
                ASSERT(this%edges(i_edge)%volumes(i_volume) > 0)
                ASSERT(this%edges(i_edge)%volumes(i_volume) <= this%nb_volumes)
            end do
!
            ASSERT(this%edges(i_edge)%nb_faces <= 20)
            do i_face = 1, this%edges(i_edge)%nb_faces
                ASSERT(this%edges(i_edge)%faces(i_face) > 0)
                ASSERT(this%edges(i_edge)%faces(i_face) <= this%nb_faces)
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
!
            ASSERT(this%faces(i_face)%nb_volumes <= 2)
            do i_volume = 1, this%faces(i_face)%nb_volumes
                ASSERT(this%faces(i_face)%volumes(i_volume) > 0)
                ASSERT(this%faces(i_face)%volumes(i_volume) <= this%nb_volumes)
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
    subroutine refine_cell(this, cell_id)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: cell_id
!
        integer(kind=8) :: cell_type, cell_dim, cell_nodes(27), nb_nodes, cell_index, cell_type_sub
        integer(kind=8) :: nb_sub, sub_type(10), sub_loc(10, 10), i_sub, i_node, obj, cell_id_sub
        integer(kind=8) :: nodes_loc(27), nno, conv_type(10), edges(12), nb_edges, i_edge, edge_id
        integer(kind=8) :: nb_sub2, sub_type2(10), sub_loc2(10, 10), conv_type2(10), i_face, face_id
!
        cell_dim = this%cells(cell_id)%dim
        cell_type = this%cells(cell_id)%type
        obj = this%cells(cell_id)%ss_id

! --- compute sub-division
        call dividing_cell(cell_type, nb_sub, sub_type, sub_loc, conv_type)
!
        this%nb_cells = this%nb_cells-1
        this%cells(cell_id)%keep = ASTER_FALSE
        if ((this%nb_total_cells+nb_sub) >= this%max_cells) then
            call this%increase_memory("CELLS   ", &
                                      max(2*this%max_cells, this%nb_total_cells+nb_sub))
        end if

        if (this%debug) then
            print *, "Refine ", cell_id, ": ", cell_type, &
                this%converter%name(cell_type), cell_dim
        end if
!
! --- we begin by refining edges
        nb_edges = 0
        if (cell_dim == 3) then
            edges = this%volumes(obj)%edges
            nb_edges = this%volumes(obj)%nb_edges
        elseif (cell_dim == 2) then
            edges(1:4) = this%faces(obj)%edges
            nb_edges = this%faces(obj)%nb_edges
        end if
!
        do i_edge = 1, nb_edges
            edge_id = edges(i_edge)
            if (this%debug) then
                print *, "Refine edge", edge_id, ": "
            end if
! --- compute sub-division
            call dividing_cell(this%edges(edge_id)%type, nb_sub2, sub_type2, sub_loc2, conv_type2)

            do i_sub = 1, nb_sub2
                cell_type_sub = sub_type2(i_sub)
                nb_nodes = this%converter%nno(cell_type_sub)
                cell_nodes = 0

                do i_node = 1, nb_nodes
                    cell_nodes(i_node) = this%edges(edge_id)%nodes(sub_loc2(i_node, i_sub))
                end do
                cell_index = this%add_edge(cell_type_sub, cell_nodes, edge_id, i_sub)
            end do
        end do
!
! --- we refine the face after
        if (cell_dim == 3) then
            do i_face = 1, this%volumes(obj)%nb_faces
                if (this%debug) then
                    print *, "Refine face", face_id, ": "
                end if
                face_id = this%volumes(obj)%faces(i_face)
! --- compute sub-division
                call dividing_cell(this%faces(face_id)%type, nb_sub2, sub_type2, sub_loc2, &
                                   conv_type2)

                do i_sub = 1, nb_sub2
                    cell_type_sub = sub_type2(i_sub)
                    nb_nodes = this%converter%nno(cell_type_sub)
                    cell_nodes = 0

                    do i_node = 1, nb_nodes
                        cell_nodes(i_node) = this%faces(face_id)%nodes(sub_loc2(i_node, i_sub))
                    end do
                    cell_index = this%add_face(cell_type_sub, cell_nodes, face_id, i_sub)
                end do
            end do
        end if
!
! --- we refine the cell
        if (this%debug) then
            print *, "Refine cell", cell_id, ": "
        end if
        do i_sub = 1, nb_sub
            cell_type_sub = sub_type(i_sub)
            nb_nodes = this%converter%nno(cell_type_sub)
            cell_nodes = 0
!
            if (cell_dim == 3) then
                do i_node = 1, nb_nodes
                    cell_nodes(i_node) = this%volumes(obj)%nodes(sub_loc(i_node, i_sub))
                end do
                cell_index = this%add_volume(cell_type_sub, cell_nodes, obj, i_sub)
            elseif (cell_dim == 2) then
                do i_node = 1, nb_nodes
                    cell_nodes(i_node) = this%faces(obj)%nodes(sub_loc(i_node, i_sub))
                end do
                cell_index = this%add_face(cell_type_sub, cell_nodes, obj, i_sub)
            elseif (cell_dim == 1) then
                do i_node = 1, nb_nodes
                    cell_nodes(i_node) = this%edges(obj)%nodes(sub_loc(i_node, i_sub))
                end do
                cell_index = this%add_edge(cell_type_sub, cell_nodes, obj, i_sub)
            else
                ASSERT(ASTER_FALSE)
            end if
!
            this%nb_total_cells = this%nb_total_cells+1
            ASSERT(this%nb_total_cells <= this%max_cells)
            this%nb_cells = this%nb_cells+1
            cell_id_sub = this%nb_total_cells
            this%cells(cell_id_sub)%dim = cell_dim
            this%cells(cell_id_sub)%id = cell_id_sub
            this%cells(cell_id_sub)%ss_id = cell_index
!
            this%cells(cell_id)%child(i_sub) = cell_id_sub
            this%cells(cell_id)%nb_child = this%cells(cell_id)%nb_child+1
!
! --- Il faut convertir comme le type de depart et pas lineaire
            this%cells(cell_id_sub)%type = conv_type(i_sub)
            nno = this%converter%nno(this%cells(cell_id_sub)%type)
            call this%numbering_nodes(this%cells(cell_id_sub)%type, nodes_loc)
            this%cells(cell_id_sub)%nodes = 0
            if (cell_dim == 3) then
                do i_node = 1, nno
                    this%cells(cell_id_sub)%nodes(i_node) = &
                        this%volumes(cell_index)%nodes(nodes_loc(i_node))
                end do
            elseif (cell_dim == 2) then
                do i_node = 1, nno
                    this%cells(cell_id_sub)%nodes(i_node) = &
                        this%faces(cell_index)%nodes(nodes_loc(i_node))
                end do
            elseif (cell_dim == 1) then
                do i_node = 1, nno
                    this%cells(cell_id_sub)%nodes(i_node) = &
                        this%edges(cell_index)%nodes(nodes_loc(i_node))
                end do
            else
                ASSERT(ASTER_FALSE)
            end if
        end do
!
    end subroutine
!
! ==================================================================================================
!
    subroutine refine(this, level)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: level
! -----------------------------------------------------------------------
        integer(kind=8) :: i_cell, nb_cells_ref, i_level
        real(kind=8) :: start, end
!
! --- Refine cells
        this%nb_level = level
        if (this%info >= 2) then
            print *, "Refining mesh..."
            call cpu_time(start)
        end if
        do i_level = 1, level
            if (this%info >= 2) then
                print *, "- Level ", i_level
            end if
            nb_cells_ref = this%nb_total_cells
            do i_cell = 1, nb_cells_ref
                if (this%cells(i_cell)%keep .and. this%cells(i_cell)%dim > 0) then
                    call this%refine_cell(i_cell)
                end if
            end do
        end do
!
        call this%update()
        if (this%info >= 2) then
            call cpu_time(end)
            print *, "... in ", end-start, " seconds."
        end if
    end subroutine
!
! ==================================================================================================
!
    recursive subroutine sub_cells(this, cell_id, nb_cells, cells)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        integer(kind=8), intent(inout) :: nb_cells, cell_id
        integer(kind=8), allocatable :: cells(:)
! -----------------------------------------------------------------------
        integer(kind=8) :: i_cell
!
        if (this%cells(cell_id)%keep) then
            nb_cells = nb_cells+1
            ASSERT(nb_cells <= size(cells))
            ASSERT(this%cells(cell_id)%nb_child == 0)
            cells(nb_cells) = this%cells(cell_id)%id
        else
            do i_cell = 1, this%cells(cell_id)%nb_child
                call this%sub_cells(this%cells(cell_id)%child(i_cell), nb_cells, cells)
            end do
        end if

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
        integer(kind=8) :: i_node, nb_nodes_in, nb_nodes_out, codret, nb_grno_out
        integer(kind=8) :: i_group, node_id, nb_grno_in
        character(len=24) :: grno_in, nomgrp
        integer(kind=8), pointer :: nodes_in(:) => null()
        integer(kind=8), pointer :: grno_out(:) => null()
        integer(kind=8), pointer :: nodes_out(:) => null()
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
                if (this%nodes(node_id)%keep) then
                    grno_out(i_group) = grno_out(i_group)+1
                end if
            end do
            if (grno_out(i_group) > 0) then
                nb_grno_out = nb_grno_out+1
            end if
        end do
!
! --- Create groups
        if (nb_grno_out > 0) then
            call jecreo(gpptnn, 'G N K24')
            call jeecra(gpptnn, 'NOMMAX', nb_grno_out)
            call jecrec(grpnoe, 'G V I', 'NO '//gpptnn, 'DISPERSE', 'VARIABLE', nb_grno_out)
!
            do i_group = 1, nb_grno_in
                if (grno_out(i_group) > 0) then
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
                        if (this%nodes(node_id)%keep) then
                            nb_nodes_out = nb_nodes_out+1
                            nodes_out(nb_nodes_out) = this%nodes(node_id)%id
                        end if
                    end do
                end if
            end do
        end if
!
        AS_DEALLOCATE(vi=grno_out)
!
999     continue
!
        call jedema()
!
    end subroutine
!
! ==================================================================================================
!
    subroutine copy_group_ma(this, grpmai, gpptnm)
!
        implicit none
!
        class(Mmesh), intent(inout) :: this
        character(len=24), intent(in) :: grpmai, gpptnm
! -----------------------------------------------------------------------
        integer(kind=8) :: i_cell, nb_cells_in, nb_cells_out, codret, nb_grma_out
        integer(kind=8) :: i_group, cell_id, nb_grma_in
        integer(kind=8) :: nb_cells
        character(len=24) :: grma_in, nomgrp
        integer(kind=8), pointer :: cells_in(:) => null()
        integer(kind=8), pointer :: grma_out(:) => null()
        integer(kind=8), pointer :: cells_out(:) => null()
        integer(kind=8), allocatable :: cells(:)
!
        call jemarq()
!
        grma_in = this%mesh_in//'.GROUPEMA'
!
        call jedetr(grpmai)
        call jedetr(gpptnm)
        call jeexin(grma_in, codret)
        if (codret .eq. 0) goto 999

        call jelira(grma_in, 'NOMUTI', nb_grma_in)
        AS_ALLOCATE(vi=grma_out, size=nb_grma_in)
        allocate (cells(8**this%nb_level))
        grma_out(:) = 0
        nb_grma_out = 0
!
! --- Find groups
        do i_group = 1, nb_grma_in
            call jeveuo(jexnum(grma_in, i_group), 'L', vi=cells_in)
            call jelira(jexnum(grma_in, i_group), 'LONUTI', nb_cells_in)
            do i_cell = 1, nb_cells_in
                cell_id = cells_in(i_cell)
                nb_cells = 0
                call this%sub_cells(cell_id, nb_cells, cells)
                grma_out(i_group) = grma_out(i_group)+nb_cells
            end do
            if (grma_out(i_group) > 0) then
                nb_grma_out = nb_grma_out+1
            end if
        end do
!
! --- Create groups
        if (nb_grma_out > 0) then
            call jecreo(gpptnm, 'G N K24')
            call jeecra(gpptnm, 'NOMMAX', nb_grma_out)
            call jecrec(grpmai, 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE', nb_grma_out)
!
            do i_group = 1, nb_grma_in
                if (grma_out(i_group) > 0) then
                    call jenuno(jexnum(grma_in, i_group), nomgrp)
                    call jecroc(jexnom(grpmai, nomgrp))
                    call jeveuo(jexnum(grma_in, i_group), 'L', vi=cells_in)
                    call jelira(jexnum(grma_in, i_group), 'LONUTI', nb_cells_in)
                    call jeecra(jexnom(grpmai, nomgrp), 'LONMAX', grma_out(i_group))
                    call jeecra(jexnom(grpmai, nomgrp), 'LONUTI', grma_out(i_group))
                    call jeveuo(jexnom(grpmai, nomgrp), 'E', vi=cells_out)
                    nb_cells_out = 1
                    do i_cell = 1, nb_cells_in
                        cell_id = cells_in(i_cell)
                        nb_cells = 0
                        call this%sub_cells(cell_id, nb_cells, cells)
                        cells_out(nb_cells_out:nb_cells_out+nb_cells) = cells(1:nb_cells)
                        nb_cells_out = nb_cells_out+nb_cells
                    end do
                end if
            end do
        end if
!
        AS_DEALLOCATE(vi=grma_out)
        deallocate (cells)
!
999     continue
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
        integer(kind=8), intent(in) :: new_size
! -----------------------------------------------------------------------
        integer(kind=8) :: old_size
        type(Mnode), allocatable :: nodes(:)
        type(Medge), allocatable :: edges(:)
        type(Mface), allocatable :: faces(:)
        type(Mvolume), allocatable :: volumes(:)
        type(Mcell), allocatable :: cells(:)
        integer(kind=4), allocatable :: lastGhostsLayer(:)
        integer(kind=8), allocatable :: list(:)
!
        if (object == "NODES") then
            old_size = size(this%nodes)
            ! nodes
            allocate (nodes(old_size))
            nodes(1:old_size) = this%nodes(1:old_size)
            deallocate (this%nodes)
            allocate (this%nodes(new_size))
            this%nodes(1:old_size) = nodes(1:old_size)
            deallocate (nodes)
            this%max_nodes = new_size
            ! ghosts
            allocate (lastGhostsLayer(old_size))
            lastGhostsLayer(1:old_size) = this%lastGhostsLayer(1:old_size)
            deallocate (this%lastGhostsLayer)
            allocate (this%lastGhostsLayer(new_size))
            this%lastGhostsLayer(1:old_size) = lastGhostsLayer(1:old_size)
            deallocate (lastGhostsLayer)
        elseif (object == "EDGES") then
            old_size = size(this%edges)
            allocate (edges(old_size))
            edges(1:old_size) = this%edges(1:old_size)
            deallocate (this%edges)
            allocate (this%edges(new_size))
            this%edges(1:old_size) = edges(1:old_size)
            deallocate (edges)
            this%max_edges = new_size
        elseif (object == "EDGES_DE") then
            old_size = size(this%edges_dege)
            allocate (list(old_size))
            list(1:old_size) = this%edges_dege(1:old_size)
            deallocate (this%edges_dege)
            allocate (this%edges_dege(new_size))
            this%edges_dege(1:old_size) = list(1:old_size)
            this%edges_dege(old_size+1:new_size) = 0
            deallocate (list)
            this%max_edges_dege = new_size
        elseif (object == "FACES") then
            old_size = size(this%faces)
            allocate (faces(old_size))
            faces(1:old_size) = this%faces(1:old_size)
            deallocate (this%faces)
            allocate (this%faces(new_size))
            this%faces(1:old_size) = faces(1:old_size)
            deallocate (faces)
            this%max_faces = new_size
        elseif (object == "FACES_DE") then
            old_size = size(this%faces_dege)
            allocate (list(old_size))
            list(1:old_size) = this%faces_dege(1:old_size)
            deallocate (this%faces_dege)
            allocate (this%faces_dege(new_size))
            this%faces_dege(1:old_size) = list(1:old_size)
            this%faces_dege(old_size+1:new_size) = 0
            deallocate (list)
            this%max_faces_dege = new_size
        elseif (object == "VOLUMES") then
            old_size = size(this%volumes)
            allocate (volumes(old_size))
            volumes(1:old_size) = this%volumes(1:old_size)
            deallocate (this%volumes)
            allocate (this%volumes(new_size))
            this%volumes(1:old_size) = volumes(1:old_size)
            deallocate (volumes)
            this%max_volumes = new_size
        elseif (object == "CELLS") then
            old_size = size(this%cells)
            allocate (cells(old_size))
            cells(1:old_size) = this%cells(1:old_size)
            deallocate (this%cells)
            allocate (this%cells(new_size))
            this%cells(1:old_size) = cells(1:old_size)
            deallocate (cells)
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
        character(len=24) :: send, recv, domj, gcom, pgin, nblg
        character(len=19) :: joints
        integer(kind=8), pointer :: v_rnode(:) => null()
        integer(kind=8), pointer :: v_noex(:) => null()
        integer(kind=8), pointer :: v_nojoin(:) => null()
        integer(kind=8), pointer :: v_snume(:) => null()
        integer(kind=8), pointer :: v_rnume(:) => null()
        integer(kind=8), pointer :: v_nkeep(:) => null()
        integer(kind=8), pointer :: v_ckeep(:) => null()
        integer(kind=8), pointer :: v_nulogl(:) => null()
        integer(kind=8), pointer :: v_proc(:) => null()
        integer(kind=4), pointer :: v_pgid(:) => null()
        integer(kind=8), pointer :: v_gcom(:) => null()
        integer(kind=8), pointer :: v_comm(:) => null()
        integer(kind=8), pointer :: v_tag(:) => null()
        integer(kind=8), pointer :: v_nblg(:) => null()
        aster_logical, pointer :: v_keep(:) => null()
        real(kind=8), pointer :: v_send(:) => null()
        real(kind=8), pointer :: v_recv(:) => null()
        mpi_int, parameter :: mpi_one = to_mpi_int(1)
        mpi_int :: msize, mrank, count_send, count_recv, id, tag, mpicou
        integer(kind=8) :: nbproc, rank, ind, nb_recv, i_proc, recv1(1)
        integer(kind=8) :: n_coor_send, n_coor_recv, proc_id, i_comm, domj_i
        integer(kind=8) :: i_node, nb_nodes_keep, i_node_r, node_id, j_node
        integer(kind=8) :: i_cell, nno, nb_cells_keep, owner, cell_id, i_layer
        real(kind=8) :: coor(3), coor_diff(3), tole_comp, start, end, coor_i(3), coor_j(3)
        real(kind=8), parameter :: tole = 1.d-15
        aster_logical :: find, keep
        ! Allows to send all nodes of a mesh to seach candidates - robust bus very slow - for debug
        aster_logical, parameter :: all_nodes = ASTER_FALSE
!
        if (isParallelMesh(mesh_out)) then
            if (this%info >= 2) then
                print *, "Create joints..."
                call cpu_time(start)
            end if

            joints = mesh_out//".JOIN"
            domj = joints//".DOMJ"
            send = joints//".SEND"
            recv = joints//".RECV"
            gcom = joints//".GCOM"
            pgin = joints//".PGID"
            nblg = joints//".NBLG"

            call jemarq()
            call asmpi_comm('GET', mpicou)
            call asmpi_info(rank=mrank, size=msize)
            rank = to_aster_int(mrank)
            nbproc = to_aster_int(msize)
!
! --- 0: On enregiste le nombre de couches de ghost
            call wkvect(nblg, 'G V I', 1, vi=v_nblg)
            v_nblg(1) = this%nb_layer
            call jedetr(pgin)
            call wkvect(pgin, 'G V S', nbproc, vi4=v_pgid)
            v_pgid(1:nbproc) = to_mpi_int(-1)
!
! --- 1: On commence par compter le nombre de noeuds que l'on doit recevoir
!
            call jeveuo(mesh_out//".NOEX", 'L', vi=v_noex)
            call jeveuo(mesh_out//".NUNOLG", 'L', vi=v_nulogl)

            if (nbproc == 1) goto 100

            call wkvect("&&CREAMA.RNODE", 'V V I', nbproc, vi=v_rnode)
            do i_node = 1, this%nb_nodes
                owner = v_noex(i_node)
                v_rnode(owner+1) = v_rnode(owner+1)+1
            end do
            v_rnode(rank+1) = 0
! --- On compte combien on doit recevoir et envoyer
!
            nb_recv = 0
            do i_proc = 0, nbproc-1
                if (v_rnode(i_proc+1) > 0) then
                    nb_recv = nb_recv+1
                end if
            end do
            if (nb_recv > 0) then
                call wkvect(domj, 'G V I', nb_recv, vi=v_proc)
                call jecrec(send, 'G V I', 'NU', 'DISPERSE', 'VARIABLE', nb_recv)
                call jecrec(recv, 'G V I', 'NU', 'DISPERSE', 'VARIABLE', nb_recv)
                call jeveuo(gcom, 'E', vi=v_gcom)
                v_gcom(1) = mpicou
            else
                goto 100
            end if
            nb_recv = 0
            do i_proc = 0, nbproc-1
                if (v_rnode(i_proc+1) > 0) then
                    nb_recv = nb_recv+1
                    v_proc(nb_recv) = i_proc
                    v_pgid(i_proc+1) = to_mpi_int(i_proc)
                else
                    v_pgid(i_proc+1) = to_mpi_int(-1)
                end if
            end do
            call sort_i8(v_proc, nb_recv)
!
            call wkvect('&&CREAMA.COMM', 'V V I', nbproc, vi=v_comm)
            call wkvect('&&CREAMA.TAG', 'V V I', nbproc, vi=v_tag)
            call build_tree_comm(v_proc, nb_recv, v_pgid, mpicou, v_comm, v_tag)
! --- Pour accélérer la recherche, on garde les cells avec un noeud non-proprio
            call wkvect('&&CREAMA.CKEEP', 'V V I', this%nb_total_cells, vi=v_ckeep)
            nb_cells_keep = 0
            do i_cell = 1, this%nb_total_cells
                if (.not. this%cells(i_cell)%keep) cycle
                nno = this%converter%nno(this%cells(i_cell)%type)
                do i_node = 1, nno
                    if (v_noex(this%nodes(this%cells(i_cell)%nodes(i_node))%id) .ne. rank) then
                        nb_cells_keep = nb_cells_keep+1
                        v_ckeep(nb_cells_keep) = i_cell
                        exit
                    end if
                end do
            end do
!
            call wkvect('&&CREAMA.NKEEP', 'V V L', this%nb_total_nodes, vl=v_keep)
            call wkvect('&&CREAMA.NRESTRICT', 'V V I', this%nb_total_nodes, vi=v_nkeep)
!
            if (this%info >= 2) then
                print *, "-Nombre de cells candidates totales: ", nb_cells_keep
            end if
! --- On crée les joints
            do i_comm = 1, nb_recv
                domj_i = v_comm(i_comm)
                proc_id = v_proc(domj_i)
                tag = to_mpi_int(v_tag(i_comm))
                id = to_mpi_int(proc_id)
!
! --- Send and Receive size
                n_coor_send = v_rnode(proc_id+1)
                call asmpi_sendrecv_i([n_coor_send], mpi_one, id, tag, &
                                      recv1, mpi_one, id, tag, mpicou)
                n_coor_recv = recv1(1)
!
                call wkvect("&&CREAMA.SCOOR", 'V V R', 4*n_coor_send, vr=v_send)
                call wkvect("&&CREAMA.RCOOR", 'V V R', 4*n_coor_recv, vr=v_recv)
!
! --- Prepare data to send (local_id, coor(1:3))
                ind = 0
                do i_node = 1, this%nb_total_nodes
                    if (this%nodes(i_node)%keep) then
                        if (v_noex(this%nodes(i_node)%id) == proc_id) then
                            v_send(ind+1) = real(this%nodes(i_node)%id, kind=8)
                            v_send(ind+2:ind+4) = this%nodes(i_node)%coor(1:3)
                            ind = ind+4
                        end if
                    end if
                end do
                ASSERT(ind == 4*n_coor_send)
!
! --- Send and Receive data
                count_send = to_mpi_int(4*n_coor_send)
                count_recv = to_mpi_int(4*n_coor_recv)
                call asmpi_sendrecv_r(v_send, count_send, id, tag, &
                                      v_recv, count_recv, id, tag, mpicou)
!
! --- Research corresponding node
! --- On cherche les noeuds avec les coor - c'est pas génial mais pas mieux
! --- Pour accélérer la recherche, on garde les noeuds voisins des non-proprio
                v_keep = ASTER_FALSE
                nb_nodes_keep = 0
!
! --- Add first layer of cells
                do i_cell = 1, nb_cells_keep
                    cell_id = v_ckeep(i_cell)
                    keep = ASTER_FALSE
                    nno = this%converter%nno(this%cells(cell_id)%type)
                    do i_node = 1, nno
                        owner = v_noex(this%nodes(this%cells(cell_id)%nodes(i_node))%id)
                        if (owner == proc_id .or. all_nodes) then
                            keep = ASTER_TRUE
                            exit
                        end if
                    end do
                    if (keep) then
                        do i_node = 1, nno
                            node_id = this%cells(cell_id)%nodes(i_node)
                            owner = v_noex(this%nodes(node_id)%id)
                            if (((owner == rank .or. all_nodes) .and. (.not. v_keep(node_id)))) then
                                v_keep(node_id) = ASTER_TRUE
                                nb_nodes_keep = nb_nodes_keep+1
                                v_nkeep(nb_nodes_keep) = node_id
                            end if
                        end do
                    end if
                end do
!
! --- Add additional layer
                do i_layer = 2, this%nb_layer
                    do i_cell = 1, this%nb_total_cells
                        if (.not. this%cells(i_cell)%keep) cycle
                        nno = this%converter%nno(this%cells(i_cell)%type)
                        keep = ASTER_FALSE
                        do i_node = 1, nno
                            node_id = this%cells(i_cell)%nodes(i_node)
                            owner = v_noex(this%nodes(node_id)%id)
                            if (owner == rank .and. v_keep(node_id)) then
                                keep = ASTER_TRUE
                                exit
                            end if
                        end do
                        if (keep) then
                            do i_node = 1, nno
                                node_id = this%cells(i_cell)%nodes(i_node)
                                owner = v_noex(this%nodes(node_id)%id)
                                if (owner == rank .and. .not. v_keep(node_id)) then
                                    v_keep(node_id) = ASTER_TRUE
                                    nb_nodes_keep = nb_nodes_keep+1
                                    v_nkeep(nb_nodes_keep) = node_id
                                end if
                            end do
                        end if
                    end do
                end do
!
                if (this%info >= 2) then
                    print *, "-Domaine: ", proc_id, &
                        ", nombre de noeuds à trouver: ", n_coor_recv, &
                        " pour ", nb_nodes_keep, " candidats"
                end if
!
                !  On regarde que les noeuds ne sont pas confondu.
                do i_node = 1, nb_nodes_keep
                    coor_i = this%nodes(v_nkeep(i_node))%coor
                    tole_comp = max(tole, tole*norm2(coor_i))
                    do j_node = i_node+1, nb_nodes_keep
                        coor_j = this%nodes(v_nkeep(j_node))%coor
                        coor_diff = abs(coor_i-coor_j)
                        if (maxval(coor_diff) < tole_comp) then
                            !! Verif pas de noeud double
                            call utmess('F', 'MAILLAGE1_4')
                        end if
                    end do
                end do

                do i_node = 1, n_coor_recv
                    coor_i = v_recv(4*(i_node-1)+2:4*(i_node-1)+4)
                    tole_comp = max(tole, tole*norm2(coor_i))
                    do j_node = i_node+1, n_coor_recv
                        coor_j = v_recv(4*(j_node-1)+2:4*(j_node-1)+4)
                        coor_diff = abs(coor_i-coor_j)
                        if (maxval(coor_diff) < tole_comp) then
                            !! Verif pas de noeud double
                            call utmess('F', 'MAILLAGE1_4')
                        end if
                    end do
                end do
!
! --- Create joint .E
                call jecroc(jexnum(send, i_comm))
                call jeecra(jexnum(send, domj_i), 'LONMAX', 2*n_coor_recv)
                call jeecra(jexnum(send, domj_i), 'LONUTI', 2*n_coor_recv)
                call jeveuo(jexnum(send, domj_i), 'E', vi=v_nojoin)
!
                call wkvect("&&CREAMA.SNUME", 'V V I', 2*n_coor_recv, vi=v_snume)
!
! --- Search nodes with coordinates
                ASSERT(nb_nodes_keep >= n_coor_recv)
                v_keep(1:nb_nodes_keep) = ASTER_TRUE
                do i_node_r = 1, n_coor_recv
                    find = ASTER_FALSE
                    coor = v_recv(4*(i_node_r-1)+2:4*(i_node_r-1)+4)
                    node_id = 0
                    tole_comp = max(tole, tole*norm2(coor))
                    do i_node = 1, nb_nodes_keep
                        if (v_keep(i_node)) then
                            coor_diff = abs(coor-this%nodes(v_nkeep(i_node))%coor)
                            if (maxval(coor_diff) < tole_comp) then
                                find = ASTER_TRUE
                                node_id = this%nodes(v_nkeep(i_node))%id
                                v_keep(i_node) = ASTER_FALSE
                                exit
                            end if
                        end if
                    end do
                    ASSERT(find)
                    v_nojoin(2*(i_node_r-1)+1) = node_id
                    v_nojoin(2*(i_node_r-1)+2) = int(v_recv(4*(i_node_r-1)+1))
                    v_snume(2*(i_node_r-1)+1) = node_id
                    v_snume(2*(i_node_r-1)+2) = v_nulogl(node_id)
                end do
!
! --- Send and recv data
                call wkvect("&&CREAMA.RNUME", 'V V I', 2*n_coor_send, vi=v_rnume)

                count_send = to_mpi_int(2*n_coor_recv)
                count_recv = to_mpi_int(2*n_coor_send)
                call asmpi_sendrecv_i(v_snume, count_send, id, tag, &
                                      v_rnume, count_recv, id, tag, mpicou)
!
! --- Create joint .R
                call jecroc(jexnum(recv, i_comm))
                call jeecra(jexnum(recv, domj_i), 'LONMAX', 2*n_coor_send)
                call jeecra(jexnum(recv, domj_i), 'LONUTI', 2*n_coor_send)
                call jeveuo(jexnum(recv, domj_i), 'E', vi=v_nojoin)
!
                do i_node = 1, n_coor_send
                    v_nojoin(2*(i_node-1)+1) = int(v_send(4*(i_node-1)+1))
                    v_nojoin(2*(i_node-1)+2) = v_rnume(2*(i_node-1)+1)
                    ASSERT(v_nulogl(v_nojoin(2*(i_node-1)+1)) == -1)
                    v_nulogl(v_nojoin(2*(i_node-1)+1)) = v_rnume(2*(i_node-1)+2)
                end do
! --- Cleaning
                call jedetr("&&CREAMA.SCOOR")
                call jedetr("&&CREAMA.SNUME")
                call jedetr("&&CREAMA.RCOOR")
                call jedetr("&&CREAMA.RNUME")
            end do
! --- Cleaning
            call jedetr('&&CREAMA.NRESTRICT')
            call jedetr('&&CREAMA.NKEEP')
            call jedetr('&&CREAMA.CKEEP')
            call jedetr("&&CREAMA.RNODE")
            call jedetr("&&CREAMA.TAG")
            call jedetr("&&CREAMA.COMM")
!
100         continue
!
            if (nbproc == 1) then
                do i_node = 1, this%nb_nodes
                    v_nulogl(i_node) = i_node
                end do
            end if
!
! --- verify
            do i_node = 1, this%nb_nodes
                ASSERT(v_nulogl(i_node) >= 0)
            end do
!
            if (this%info >= 2) then
                call cpu_time(end)
                print *, "... in ", end-start, " seconds."
            end if
!
            call jedema()
        end if
    end subroutine
!
! ==================================================================================================
!
    logical function mult_elem(nb_nodes, nodes)
!
! Performs circular permutation such that the first element is the smallest and
! the second element is the second smallest
!
        implicit none
!
        integer(kind=8), intent(in) :: nb_nodes
        integer(kind=8), intent(inout) :: nodes(1:nb_nodes)
!
        integer(kind=8) :: i_node, j_node
!
        ASSERT(nb_nodes <= 27)
        mult_elem = .false.
!
        do i_node = 1, nb_nodes
            do j_node = i_node+1, nb_nodes
                if (nodes(i_node) == nodes(j_node)) then
                    mult_elem = .true.
                    exit
                end if
            end do
        end do
    end function
!
! ==================================================================================================
!
    subroutine fix_measure(this, cell_id, cell_dim)
!
        implicit none
!
        type(Mmesh), intent(inout) :: this
        integer(kind=8), intent(in) :: cell_id, cell_dim
!
        integer(kind=8) :: i, nodes(27), nno, nbpg, ndim
        real(kind=8) :: jaco(3, 3), jacob, coopg(3), dbasis(3, 8), coorno(3), poipg(2)
        character(len=3) :: typema
! ---------------------------------------------------------------------------------
!
        if (cell_dim == 3) then
            select case (this%volumes(cell_id)%type)
            case (MT_HEXA27)
                typema = "HE8"
            case (MT_TETRA15)
                typema = "TE4"
            case (MT_PENTA21)
                typema = "PY5"
            case (MT_PYRAM19)
                typema = "PE6"
            case default
                ASSERT(ASTER_FALSE)
            end select
!
            dbasis = 0.d0
            call elraga(typema, "FPG1", ndim, nbpg, coopg, poipg)
            ASSERT(nbpg == 1)
            call elrfdf(typema, coopg, dbasis, nno_=nno)
!
! ---  Compute the jacobienne
            jaco = 0.d0
            do i = 1, nno
                coorno = this%nodes(this%volumes(cell_id)%nodes(i))%coor
                jaco(1:3, 1) = jaco(1:3, 1)+coorno(1)*dbasis(1:3, i)
                jaco(1:3, 2) = jaco(1:3, 2)+coorno(2)*dbasis(1:3, i)
                jaco(1:3, 3) = jaco(1:3, 3)+coorno(3)*dbasis(1:3, i)
            end do
!
            jacob = jaco(1, 1)*jaco(2, 2)*jaco(3, 3)+jaco(1, 3)*jaco(2, 1)*jaco(3, 2) &
                    +jaco(3, 1)*jaco(1, 2)*jaco(2, 3)-jaco(3, 1)*jaco(2, 2)*jaco(1, 3) &
                    -jaco(3, 3)*jaco(2, 1)*jaco(1, 2)-jaco(1, 1)*jaco(2, 3)*jaco(3, 2)
!
            if (jacob < 0.d0) then
                nodes = this%volumes(cell_id)%nodes
                select case (this%volumes(cell_id)%type)
                case (MT_HEXA27)
                    this%volumes(cell_id)%nodes(1:4) = [nodes(1), nodes(4), nodes(3), nodes(2)]
                    this%volumes(cell_id)%nodes(5:8) = [nodes(5), nodes(8), nodes(7), nodes(6)]
                    this%volumes(cell_id)%nodes(9:12) = [nodes(12), nodes(11), nodes(10), nodes(9)]
                    this%volumes(cell_id)%nodes(13:16) = [nodes(13), nodes(16), nodes(15), &
                                                          nodes(14)]
                    this%volumes(cell_id)%nodes(17:20) = [nodes(20), nodes(19), nodes(18), &
                                                          nodes(17)]
                    this%volumes(cell_id)%nodes(21:26) = [nodes(21), nodes(25), nodes(24), &
                                                          nodes(23), nodes(22), nodes(26)]
                    this%volumes(cell_id)%nodes(27) = nodes(27)
                case (MT_TETRA15)
                    this%volumes(cell_id)%nodes(1:4) = [nodes(1), nodes(3), nodes(2), nodes(4)]
                    this%volumes(cell_id)%nodes(5:7) = [nodes(7), nodes(6), nodes(5)]
                    this%volumes(cell_id)%nodes(8:10) = [nodes(8), nodes(9), nodes(10)]
                    this%volumes(cell_id)%nodes(11:15) = nodes(11:15)
                case (MT_PENTA21)
                    this%volumes(cell_id)%nodes(1:3) = [nodes(1), nodes(3), nodes(2)]
                    this%volumes(cell_id)%nodes(4:6) = [nodes(4), nodes(6), nodes(5)]
                    this%volumes(cell_id)%nodes(7:9) = [nodes(9), nodes(8), nodes(7)]
                    this%volumes(cell_id)%nodes(10:12) = [nodes(10), nodes(12), nodes(11)]
                    this%volumes(cell_id)%nodes(13:15) = [nodes(15), nodes(14), nodes(13)]
                    this%volumes(cell_id)%nodes(16:18) = [nodes(18), nodes(17), nodes(16)]
                    this%volumes(cell_id)%nodes(19:21) = nodes(19:21)
                case (MT_PYRAM19)
                    this%volumes(cell_id)%nodes(1:4) = [nodes(1), nodes(4), nodes(3), nodes(2)]
                    this%volumes(cell_id)%nodes(5) = nodes(5)
                    this%volumes(cell_id)%nodes(6:9) = [nodes(9), nodes(8), nodes(7), nodes(6)]
                    this%volumes(cell_id)%nodes(7:9) = [nodes(9), nodes(8), nodes(7)]
                    this%volumes(cell_id)%nodes(10:13) = [nodes(10), nodes(13), nodes(12), &
                                                          nodes(11)]
                    this%volumes(cell_id)%nodes(14) = nodes(14)

                    this%volumes(cell_id)%nodes(15:18) = [nodes(18), nodes(17), nodes(16), &
                                                          nodes(15)]
                    this%volumes(cell_id)%nodes(19) = nodes(19)
                case default
                    ASSERT(ASTER_FALSE)
                end select
!
                nno = this%converter%nno(this%cells(this%volumes(cell_id)%cell_id)%type)
                this%cells(this%volumes(cell_id)%cell_id)%nodes(1:nno) = &
                    this%volumes(cell_id)%nodes(1:nno)
!
! ------ Verification that volume is positive
                jaco = 0.d0
                do i = 1, nno
                    coorno = this%nodes(this%volumes(cell_id)%nodes(i))%coor
                    jaco(1:3, 1) = jaco(1:3, 1)+coorno(1)*dbasis(1:3, i)
                    jaco(1:3, 2) = jaco(1:3, 2)+coorno(2)*dbasis(1:3, i)
                    jaco(1:3, 3) = jaco(1:3, 3)+coorno(3)*dbasis(1:3, i)
                end do
                !
                jacob = jaco(1, 1)*jaco(2, 2)*jaco(3, 3)+jaco(1, 3)*jaco(2, 1)*jaco(3, 2) &
                        +jaco(3, 1)*jaco(1, 2)*jaco(2, 3)-jaco(3, 1)*jaco(2, 2)*jaco(1, 3) &
                        -jaco(3, 3)*jaco(2, 1)*jaco(1, 2)-jaco(1, 1)*jaco(2, 3)*jaco(3, 2)
                ASSERT(jacob > 0.d0)
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end subroutine
!
end module
