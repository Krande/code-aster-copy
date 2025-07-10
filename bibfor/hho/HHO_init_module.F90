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
module HHO_init_module
!
    use HHO_type
    use HHO_geometry_module
    use HHO_measure_module
    use HHO_quadrature_module
    use HHO_basis_module
    use HHO_utils_module
    use HHO_inertia_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/rcvala.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! generic routines to initialize data for HHO model
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: hhoInfoInitCell, hhoInfoInitFace
    private :: hhoGeomData, hhoGeomFace, hhoDataInit, hhoFaceInit, hhoCellInit
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGeomData(nodes_coor, nbnodes, typma, elem_dim)
!
        implicit none
!
        real(kind=8), dimension(3, 27), intent(out)  :: nodes_coor
        integer(kind=8), intent(out)                        :: nbnodes
        character(len=8), intent(out)               :: typma
        integer(kind=8), intent(out)                        :: elem_dim
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Get geometric data of the current element
!
! --------------------------------------------------------------------------------------------------
!
! Out nodes_coor        : coordinates of the nodes
! Out nbnodes           : number of nodes
! Out typma             : type of the element
! Out elem_dim          : dimension of the element
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        integer(kind=8) :: jv_geom
        integer(kind=8) :: inode, idim, iret
! --------------------------------------------------------------------------------------------------
!
! --- Init
        nodes_coor = 0.d0
        nbnodes = 0
        typma = ''
        elem_dim = 0
!
        call teattr('S', 'TYPMA', typma, iret)
        ASSERT(iret == 0)
!
        if (typma == 'H27') then
            typma = 'HEXA27'
            nbnodes = 27
            elem_dim = 3
        elseif (typma == 'QU8') then
            typma = 'QUAD8'
            nbnodes = 8
            elem_dim = 2
        elseif (typma == 'QU9') then
            typma = 'QUAD9'
            nbnodes = 9
            elem_dim = 2
        elseif (typma == 'TR6') then
            typma = 'TRIA6'
            nbnodes = 6
            elem_dim = 2
        elseif (typma == 'TR7') then
            typma = 'TRIA7'
            nbnodes = 7
            elem_dim = 2
        elseif (typma == 'T15') then
            typma = 'TETRA15'
            nbnodes = 15
            elem_dim = 3
        elseif (typma == 'P19') then
            typma = 'PYRAM19'
            nbnodes = 19
            elem_dim = 3
        elseif (typma == 'P21') then
            typma = 'PENTA21'
            nbnodes = 21
            elem_dim = 3
        else
            ASSERT(ASTER_FALSE)
        end if
!
        ASSERT((nbnodes .ge. 4) .and. (nbnodes .le. 27))
        ASSERT((elem_dim .eq. 2) .or. (elem_dim .eq. 3))
!
        call jevech('PGEOMER', 'L', jv_geom)
!
! --- Get coordinates
!
        do inode = 1, nbnodes
            do idim = 1, elem_dim
                nodes_coor(idim, inode) = zr(jv_geom+(inode-1)*elem_dim+idim-1)
            end do
        end do
!
        if (l_debug) then
            write (6, *) "hhoGeomData debug"
            write (6, *) "typma: ", typma
            write (6, *) "nbnodes: ", nbnodes
            do inode = 1, nbnodes
                write (6, *) "node ", inode, ": ", nodes_coor(1:3, inode)
            end do
            write (6, *) "end hhoGeomData debug"
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGeomFace(nodes_coor, nbnodes, typma, elem_dim)
!
        implicit none
!
        real(kind=8), dimension(3, 9), intent(out)  :: nodes_coor
        integer(kind=8), intent(out)                        :: nbnodes
        character(len=8), intent(out)               :: typma
        integer(kind=8), intent(out)                        :: elem_dim
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Get geometric data of the current element
!
! --------------------------------------------------------------------------------------------------
!
! Out nodes_coor        : coordinates of the nodes
! Out nbnodes           : number of nodes
! Out typma             : type of the element
! Out elem_dim          : dimension of the element
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        integer(kind=8) :: jv_geom
        integer(kind=8) :: inode, idim, iret
!
! --- Init
        nodes_coor = 0.d0
        nbnodes = 0
        typma = ''
        elem_dim = 0
!
        call teattr('S', 'TYPMA', typma, iret)
        ASSERT(iret == 0)
!
        if (typma == 'QU9') then
            typma = 'QUAD4'
            nbnodes = 4
            elem_dim = 2
        elseif (typma == 'TR7') then
            typma = 'TRIA3'
            nbnodes = 3
            elem_dim = 2
        elseif (typma == 'SE3') then
            typma = 'SEG2'
            nbnodes = 2
            elem_dim = 1
        else
            ASSERT(ASTER_FALSE)
        end if
!
        ASSERT((nbnodes .ge. 2) .and. (nbnodes .le. 9))
        ASSERT((elem_dim .eq. 1) .or. (elem_dim .eq. 2))
!
        call jevech('PGEOMER', 'L', jv_geom)
!
! - Get coordinates
!
        do inode = 1, nbnodes
            do idim = 1, elem_dim+1
                nodes_coor(idim, inode) = zr(jv_geom+(inode-1)*(elem_dim+1)+idim-1)
            end do
        end do
!
        if (l_debug) then
            write (6, *) "hhoGeomFace debug"
            write (6, *) "typma: ", typma
            write (6, *) "nbnodes: ", nbnodes
            do inode = 1, nbnodes
                write (6, *) "node ", inode, ": ", nodes_coor(1:3, inode)
            end do
            write (6, *) "end hhoGeomFace debug"
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoFaceInit(hhoFace, typma, ndim, nbnodes, nodes_coor, &
                           num_nodes_loc, num_face_loc, barycenter_cell)
!
        implicit none
!
        integer(kind=8), intent(in)                             :: typma
        integer(kind=8), intent(in)                             :: ndim
        real(kind=8), dimension(3, 4), intent(in)       :: nodes_coor
        integer(kind=8), intent(in)                             :: nbnodes
        integer(kind=8), dimension(5), intent(in)               :: num_nodes_loc
        integer(kind=8), intent(in)                             :: num_face_loc
        real(kind=8), dimension(3), optional, intent(in):: barycenter_cell
        type(HHO_Face), intent(out)                     :: hhoFace
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Initialize a HHO Face
!
! --------------------------------------------------------------------------------------------------
!
! In typma              : type of element
! In ndim               : dimension of the element
! In nodes_coor         : coordinates of nodes
! In nbnodes            : number of nodes
! In barycenter_cell    : barycenter of the cell
! Out hhoFace           : a HHO Face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: numsorted(4), ino
        ASSERT(nbnodes .le. 4)
!
! --- Init: be carefull the order to call the functions is important
!
        hhoFace%typema = typma
        hhoFace%nbnodes = nbnodes
        hhoFace%ndim = ndim
        hhoFace%coorno = hhoFaceInitCoor(nodes_coor, nbnodes, ndim, numsorted)
        hhoFace%barycenter = barycenter(hhoFace%coorno, hhoFace%nbnodes)
        if (present(barycenter_cell)) then
            hhoFace%normal = hhoNormalFace3(hhoFace, barycenter_cell)
        else
            hhoFace%normal = hhoNormalFace2(typma, nodes_coor)
        end if
        hhoFace%measure = hhoMeasureFace(hhoFace)
        hhoFace%diameter = hhoDiameterFace(hhoFace)
        hhoFace%axes = hhoLocalAxesFace(hhoFace)
        hhoFace%face_loc = num_face_loc
        do ino = 1, hhoFace%nbnodes
            hhoFace%nodes_loc(ino) = num_nodes_loc(numsorted(ino))
        end do
        ! Last node is the index of barycenter
        hhoFace%node_bar_loc = num_nodes_loc(hhoFace%nbnodes+1)
        hhoFace%l_jaco_cst = hhoIsJacobCst(hhoFace%typema, hhoFace%coorno, hhoFace%ndim+1)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCellInit(typma, nodes_coor, hhoCell)
!
        implicit none
!
        character(len=8), intent(in)                :: typma
        real(kind=8), dimension(3, 27), intent(in)   :: nodes_coor
        type(HHO_Cell), intent(out)                 :: hhoCell
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tools
!
! Initialize a HHO cell
!
! --------------------------------------------------------------------------------------------------
!
! In typma              : type of element
! In nodes_coor         : coordinates of nodes
! Out hhoCell           : a HHO cell
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        integer(kind=8), parameter                      :: max_faces = 6
        integer(kind=8), parameter                      :: max_nodes = 5
        integer(kind=8), dimension(max_nodes, max_faces) :: nodes_faces
        integer(kind=8), dimension(max_faces)           :: nbnodes_faces
        integer(kind=8), dimension(max_faces)           :: type_faces
        real(kind=8), dimension(3, max_nodes)    :: coor_face
        integer(kind=8) :: i_face, i_node
        aster_logical :: l_jaco
! --------------------------------------------------------------------------------------------------
! --- Init
        nodes_faces = 0
        nbnodes_faces = 0
        type_faces = 0
        coor_face = 0.d0
!
        if (typma == 'HEXA27') then
            hhoCell%typema = MT_HEXA8
            hhoCell%nbnodes = 8
            hhoCell%ndim = 3
            hhoCell%nbfaces = 6
            hhoCell%node_bar_loc = 27
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N21
            nodes_faces(1:5, 1) = (/1, 4, 3, 2, 21/)
            nbnodes_faces(1) = 4
            type_faces(1) = MT_QUAD4
! ----- Face 2 -> N22
            nodes_faces(1:5, 2) = (/1, 2, 6, 5, 22/)
            nbnodes_faces(2) = 4
            type_faces(2) = MT_QUAD4
! ----- Face 3 -> N23
            nodes_faces(1:5, 3) = (/2, 3, 7, 6, 23/)
            nbnodes_faces(3) = 4
            type_faces(3) = MT_QUAD4
! ----- Face 4 -> N24
            nodes_faces(1:5, 4) = (/3, 4, 8, 7, 24/)
            nbnodes_faces(4) = 4
            type_faces(4) = MT_QUAD4
! ----- Face 5 -> N25
            nodes_faces(1:5, 5) = (/1, 5, 8, 4, 25/)
            nbnodes_faces(5) = 4
            type_faces(5) = MT_QUAD4
! ----- Face 6 -> N26
            nodes_faces(1:5, 6) = (/5, 6, 7, 8, 26/)
            nbnodes_faces(6) = 4
            type_faces(6) = MT_QUAD4
!
        else if (typma == 'TETRA15') then
            hhoCell%typema = MT_TETRA4
            hhoCell%nbnodes = 4
            hhoCell%ndim = 3
            hhoCell%nbfaces = 4
            hhoCell%node_bar_loc = 15
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N11
            nodes_faces(1:4, 1) = (/1, 3, 2, 11/)
            nbnodes_faces(1) = 3
            type_faces(1) = MT_TRIA3
! ----- Face 2 -> N12
            nodes_faces(1:4, 2) = (/1, 2, 4, 12/)
            nbnodes_faces(2) = 3
            type_faces(2) = MT_TRIA3
! ----- Face 3 -> N13
            nodes_faces(1:4, 3) = (/1, 4, 3, 13/)
            nbnodes_faces(3) = 3
            type_faces(3) = MT_TRIA3
! ----- Face 4 -> N14
            nodes_faces(1:4, 4) = (/2, 3, 4, 14/)
            nbnodes_faces(4) = 3
            type_faces(4) = MT_TRIA3
!
        else if (typma == 'PYRAM19') then
            hhoCell%typema = MT_PYRAM5
            hhoCell%nbnodes = 5
            hhoCell%ndim = 3
            hhoCell%nbfaces = 5
            hhoCell%node_bar_loc = 19
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N14
            nodes_faces(1:5, 1) = (/1, 4, 3, 2, 14/)
            nbnodes_faces(1) = 4
            type_faces(1) = MT_QUAD4
! ----- Face 2 -> N15
            nodes_faces(1:4, 2) = (/1, 2, 5, 15/)
            nbnodes_faces(2) = 3
            type_faces(2) = MT_TRIA3
! ----- Face 3 -> N16
            nodes_faces(1:4, 3) = (/2, 3, 5, 16/)
            nbnodes_faces(3) = 3
            type_faces(3) = MT_TRIA3
! ----- Face 4 -> N17
            nodes_faces(1:4, 4) = (/3, 4, 5, 17/)
            nbnodes_faces(4) = 3
            type_faces(4) = MT_TRIA3
! ----- Face 5 -> N18
            nodes_faces(1:4, 5) = (/4, 1, 5, 18/)
            nbnodes_faces(5) = 3
            type_faces(5) = MT_TRIA3
!
        else if (typma == 'PENTA21') then
            hhoCell%typema = MT_PENTA6
            hhoCell%nbnodes = 6
            hhoCell%ndim = 3
            hhoCell%nbfaces = 5
            hhoCell%node_bar_loc = 21
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N16
            nodes_faces(1:5, 1) = (/1, 2, 5, 4, 16/)
            nbnodes_faces(1) = 4
            type_faces(1) = MT_QUAD4
! ----- Face 2 -> N17
            nodes_faces(1:5, 2) = (/2, 3, 6, 5, 17/)
            nbnodes_faces(2) = 4
            type_faces(2) = MT_QUAD4
! ----- Face 3 -> N18
            nodes_faces(1:5, 3) = (/1, 4, 6, 3, 18/)
            nbnodes_faces(3) = 4
            type_faces(3) = MT_QUAD4
! ----- Face 4 -> N19
            nodes_faces(1:4, 4) = (/1, 3, 2, 19/)
            nbnodes_faces(4) = 3
            type_faces(4) = MT_TRIA3
! ----- Face 5 -> N20
            nodes_faces(1:4, 5) = (/4, 5, 6, 20/)
            nbnodes_faces(5) = 3
            type_faces(5) = MT_TRIA3
!
        else if (typma == 'QUAD9') then
            hhoCell%typema = MT_QUAD4
            hhoCell%nbnodes = 4
            hhoCell%ndim = 2
            hhoCell%nbfaces = 4
            hhoCell%node_bar_loc = 9
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N5
            nodes_faces(1:3, 1) = (/1, 2, 5/)
            nbnodes_faces(1) = 2
            type_faces(1) = MT_SEG2
! ----- Face 2 -> N6
            nodes_faces(1:3, 2) = (/2, 3, 6/)
            nbnodes_faces(2) = 2
            type_faces(2) = MT_SEG2
! ----- Face 3 -> N7
            nodes_faces(1:3, 3) = (/3, 4, 7/)
            nbnodes_faces(3) = 2
            type_faces(3) = MT_SEG2
! ----- Face 4 -> N8
            nodes_faces(1:3, 4) = (/4, 1, 8/)
            nbnodes_faces(4) = 2
            type_faces(4) = MT_SEG2
!
        else if (typma == 'TRIA7') then
            hhoCell%typema = MT_TRIA3
            hhoCell%nbnodes = 3
            hhoCell%ndim = 2
            hhoCell%nbfaces = 3
            hhoCell%node_bar_loc = 7
!
! ----- !!!! Attention l'ordre des faces doit etre le meme que celui du catalogue
        !!!! Le dernier noeud est le barycentre
! ----- Face 1 -> N4
            nodes_faces(1:3, 1) = (/1, 2, 4/)
            nbnodes_faces(1) = 2
            type_faces(1) = MT_SEG2
! ----- Face 2 -> N5
            nodes_faces(1:3, 2) = (/2, 3, 5/)
            nbnodes_faces(2) = 2
            type_faces(2) = MT_SEG2
! ----- Face 3 -> N6
            nodes_faces(1:3, 3) = (/3, 1, 6/)
            nbnodes_faces(3) = 2
            type_faces(3) = MT_SEG2
!
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ----- Copy coordinates
        hhoCell%coorno = nodes_coor
        hhoCell%barycenter = barycenter(hhoCell%coorno, hhoCell%nbnodes)
        hhoCell%measure = hhoMeasureCell(hhoCell)
        hhoCell%diameter = hhoDiameterCell(hhoCell)
!
        hhoCell%axes = hhoLocalAxesCell(hhoCell)
!
! ----- Init faces
!
        do i_face = 1, hhoCell%nbfaces
            coor_face = 0.d0
            do i_node = 1, nbnodes_faces(i_face)
                coor_face(1:3, i_node) = hhoCell%coorno(1:3, nodes_faces(i_node, i_face))
            end do
            call hhoFaceInit(hhoCell%faces(i_face), type_faces(i_face), hhoCell%ndim-1, &
                             nbnodes_faces(i_face), coor_face, nodes_faces(:, i_face), &
                             i_face, hhoCell%barycenter)
        end do
!
! ----- Jacobienne
!
        if (hhoCell%typema == MT_TETRA4 .or. hhoCell%typema == MT_TRIA3) then
            hhoCell%l_jaco_cst = ASTER_TRUE
        else
            hhoCell%l_jaco_cst = ASTER_TRUE
            l_jaco = ASTER_TRUE
            do i_face = 1, hhoCell%nbfaces
                if (.not. hhoCell%faces(i_face)%l_jaco_cst) then
                    hhoCell%l_jaco_cst = ASTER_FALSE
                    l_jaco = ASTER_FALSE
                    exit
                end if
            end do
            if (l_jaco) then
                hhoCell%l_jaco_cst = hhoIsJacobCst(hhoCell%typema, hhoCell%coorno, hhoCell%ndim)
            end if
        end if
!
        if (l_debug) then
            call hhoCell%print()
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoDataInit(hhoData)
!
        implicit none
!
        type(HHO_Data), intent(out)            :: hhoData
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Get hho data on the modelisation
!
! --------------------------------------------------------------------------------------------------
!
! Out hhoData  :: information on the modelisation
! --------------------------------------------------------------------------------------------------
!
        aster_logical, parameter :: l_debug = ASTER_FALSE
        aster_logical :: l_adapt_coeff
        integer(kind=8) :: grad_deg, face_deg, cell_deg, jv_mater
        integer(kind=8) ::  jtab(1), iret, icodre(1)
        real(kind=8) :: coef_stab
        real(kind=8) :: resu_vale(1)
        character(len=16), parameter :: resu_name(1) = (/"COEF_STAB"/)
!
! --- Read Parameters
!
        ASSERT(lteatt('TYPMOD2', 'HHO') .or. lteatt('TYPMOD2', 'HHO_GRAD'))
        if (lteatt('FORMULATION', 'HHO_CSTE')) then
            face_deg = 0
            cell_deg = 0
            grad_deg = 0
        elseif (lteatt('FORMULATION', 'HHO_LINE')) then
            face_deg = 1
            cell_deg = 1
            grad_deg = 1
        elseif (lteatt('FORMULATION', 'HHO_QUAD')) then
            face_deg = 2
            cell_deg = 2
            grad_deg = 2
        elseif (lteatt('FORMULATION', 'HHO_CUBI')) then
            face_deg = 3
            cell_deg = 3
            grad_deg = 3
        elseif (lteatt('FORMULATION', 'HHO_QUAR')) then
            face_deg = 4
            cell_deg = 4
            grad_deg = 4
        elseif (lteatt('FORMULATION', 'HHO_QUIN')) then
            face_deg = 5
            cell_deg = 5
            grad_deg = 5
        elseif (lteatt('FORMULATION', 'HHO_MCSTE')) then
            face_deg = 0
            cell_deg = 1
            grad_deg = 0
        elseif (lteatt('FORMULATION', 'HHO_MLINE')) then
            face_deg = 1
            cell_deg = 2
            grad_deg = 1
        elseif (lteatt('FORMULATION', 'HHO_MQUAD')) then
            face_deg = 2
            cell_deg = 3
            grad_deg = 2
        elseif (lteatt('FORMULATION', 'HHO_MCUBI')) then
            face_deg = 3
            cell_deg = 4
            grad_deg = 3
        elseif (lteatt('FORMULATION', 'HHO_MQUAR')) then
            face_deg = 4
            cell_deg = 5
            grad_deg = 4
        elseif (lteatt('FORMULATION', 'HHO_MQUIN')) then
            face_deg = 5
            cell_deg = 6
            grad_deg = 5
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Get coefficient (if possible)
!
        coef_stab = 0.d0
        iret = -1
        l_adapt_coeff = ASTER_TRUE
        call tecach('NNN', 'PMATERC', 'L', iret, nval=1, itab=jtab)
        if (iret .eq. 0) then
            call jevech('PMATERC', 'L', jv_mater)
            call rcvala(zi(jv_mater), ' ', 'HHO', 0, ' ', [0.0], &
                        1, resu_name, resu_vale, icodre, 0)
            l_adapt_coeff = icodre(1) == 1
            if (.not. l_adapt_coeff) then
                coef_stab = resu_vale(1)
            end if
        end if
!
! --- Init
!
        call hhoData%initialize(face_deg, cell_deg, grad_deg, coef_stab, l_debug, l_adapt_coeff)
!
        if (hhoData%debug()) then
            call hhoData%print()
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuad)
!
        implicit none
!
        type(HHO_Cell), intent(out)                 :: hhoCell
        type(HHO_Data), intent(out)                 :: hhoData
        integer(kind=8), intent(in), optional               :: npg
        type(HHO_Quadrature), optional, intent(out) :: hhoQuad

!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Get hho data on the modelisation
!
! --------------------------------------------------------------------------------------------------
!
!   Initialize information for a HHO cell
!   Out hhoCell  : the current HHO Cell
!   Out hhoData  : information on HHO methods
!   Out hhoQuad  : Quadrature for the cell
!   In npg (optional)      : number of quadrature point for the face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: nbnodes, elem_dim
        character(len=8) :: typma
        real(kind=8) :: coor(3, 27)
        aster_logical :: laxis
!
        coor = 0.d0
!
! --- Get HHO informations
!
        call hhoGeomData(coor, nbnodes, typma, elem_dim)
        call hhoDataInit(hhoData)
!
! --- Initialize HHO Cell
        call hhoCellInit(typma, coor, hhoCell)
!
! --- Get quadrature (optional)
!
        if (present(hhoQuad)) then
            laxis = lteatt("TYPMOD", "AXIS")
            if (present(npg)) then
                call hhoQuad%initCell(hhoCell, npg, laxis)
            else
                call hhoQuad%getQuadCell(hhoCell, 2*hhoData%cell_degree(), &
                                         laxis, param=ASTER_TRUE)
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoInfoInitFace(hhoFace, hhoData, npg, hhoQuadFace)
!
        implicit none
!
        type(HHO_Face), intent(out)   :: hhoFace
        type(HHO_Data), intent(out)   :: hhoData
        integer(kind=8), intent(in), optional :: npg
        type(HHO_Quadrature), intent(out), optional :: hhoQuadFace
!
! --------------------------------------------------------------------------------------------------
!
! HHO - generic tool
!
! Get hho info on the modelisation and initialize
!
! --------------------------------------------------------------------------------------------------
!
!   Initialize information for surfacique load
!   Out hhoFace  : the current HHO Face
!   Out hhoData  : information on HHO methods
!   Out hhoQuad  : Quadrature for the face
!   In npg (optional)      : number of quadrature point for the face
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: nbnodes, elem_dim, numnodes(5), enumf
        real(kind=8) :: nodes_coor(3, 9)
        character(len=8) :: typma
        aster_logical :: laxis
!
! --- Get HHO informations
!
        call hhoGeomFace(nodes_coor, nbnodes, typma, elem_dim)
        ASSERT(elem_dim == 1 .or. elem_dim == 2)
        call hhoDataInit(hhoData)
!
        if (typma == 'QUAD4') then
            numnodes(1:5) = (/1, 2, 3, 4, 9/)
            enumf = MT_QUAD4
        else if (typma == 'TRIA3') then
            numnodes(1:4) = (/1, 2, 3, 7/)
            enumf = MT_TRIA3
        else if (typma == 'SEG2 ') then
            numnodes(1:3) = (/1, 2, 3/)
            enumf = MT_SEG2
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Initialize HHO Face
!
        call hhoFaceInit(hhoFace, enumf, elem_dim, nbnodes, nodes_coor, numnodes, 1)

! --- Get quadrature (optional)
!
        if (present(hhoQuadFace)) then
            laxis = lteatt("TYPMOD", "AXIS")
            if (present(npg)) then
                call hhoQuadFace%initFace(hhoFace, npg, laxis)
            else
                call hhoQuadFace%GetQuadFace(hhoFace, 2*hhoData%face_degree(), &
                                             laxis, param=ASTER_TRUE)
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
