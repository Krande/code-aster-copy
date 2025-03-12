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
module mesh_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    private :: getSignNormalSkinToSupport
    public :: getPropertiesOfListOfCells, getPropertiesOfCell, getSkinCellSupport
    public :: checkNormalOnSkinCell, checkInclude, getCellOptionForName, createNameOfCell
    public :: getNodeOptionForName, createNameOfNode, getGroupsFromCell, getMeshDimension
    public :: getFirstNodeFromNodeGroup, getListOfCellGroup, checkCellsAreSkin
! ==================================================================================================
    private
#include "asterf_types.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmasu.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "blas/ddot.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/char8_to_int.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getPropertiesOfListOfCells
!
! Get properties of list of cells
!
! In  mesh             : name of mesh
! In  nbCell           : number of skin cells
! In  listCell         : list of skin cells (number)
! IO  cellNbNode       : number of nodes for each cell
! IO  cellNodeIndx     : index of first node for each cell
! Out lCellSurf        : 2D cells (surfacic) exist
! Out lCellLine        : 1D cells (lineic) exist
!
! --------------------------------------------------------------------------------------------------
    subroutine getPropertiesOfListOfCells(meshz, nbCell, listCell, cellNbNode, cellNodeIndx, &
                                          lCellSurf, lCellLine)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshz
        integer, intent(in) :: nbCell, listCell(*)
        integer, intent(inout) :: cellNbNode(nbCell), cellNodeIndx(nbCell)
        aster_logical, intent(out) :: lCellSurf, lCellLine
! - Local
        integer :: iCell, cellNume, cellTypeNume
        character(len=8) :: mesh, cellTypeName
        character(len=4) :: cellTopo
        integer, pointer :: skinCellType(:) => null()
        integer, pointer :: meshConnexLen(:) => null()
!   ------------------------------------------------------------------------------------------------
        mesh = meshz
        lCellSurf = ASTER_FALSE
        lCellLine = ASTER_FALSE
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=skinCellType)
        call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=meshConnexLen)
        do iCell = 1, nbCell
            cellNume = listCell(iCell)
            cellNbNode(iCell) = meshConnexLen(cellNume+1)-meshConnexLen(cellNume)
            cellNodeIndx(iCell) = meshConnexLen(cellNume)
            cellTypeNume = skinCellType(cellNume)
            call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
            call getPropertiesOfCell(cellTypeName, cellTopo)
            if (cellTopo .eq. 'SURF') then
                lCellSurf = ASTER_TRUE
            else if (cellTopo .eq. 'LINE') then
                lCellLine = ASTER_TRUE
            else
! --------- Ignore
            end if
        end do
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getSkinCellSupport
!
! Get "volumic" cells support of skin cells
!
! In  mesh             : name of mesh
! In  nbSkinCell       : number of skin cells
! Ptr cellSkinNume     : list of skin cells (number)
! In  lCell2d          : flag if 2d cells exist
! In  lCell1d          : flag if 1D cells exist
! Ptr cellSuppNume  : volumic" cells support of skin cells
!
! --------------------------------------------------------------------------------------------------
    subroutine getSkinCellSupport(meshz, nbSkinCell, cellSkinNume, lCell2d, lCell1d, &
                                  cellSuppNume, nbCellSupport_, suppNume_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshz
        integer, intent(in) :: nbSkinCell
        integer, pointer :: cellSkinNume(:)
        aster_logical, intent(in) :: lCell2d, lCell1d
        integer, pointer :: cellSuppNume(:)
        integer, optional, intent(in) :: nbCellSupport_
        integer, optional, pointer :: suppNume_(:)
! - Local
        character(len=8) :: mesh
        character(len=2) :: kdim
        real(kind=8), pointer :: meshNodeCoor(:) => null()
        integer, parameter :: zero = 0
        aster_logical, parameter :: skinInsideVolume = ASTER_TRUE
        integer :: ibid(1), iSkinCell
        character(len=24), parameter :: jvCellVolume = '&&UTMASU'
        integer, pointer :: vCellVolume(:) => null()
!   ------------------------------------------------------------------------------------------------
        mesh = meshz
        kdim = '  '
        if (lCell1d) then
            kdim = '2D'
            ASSERT(.not. lCell2d)
        end if
        if (lCell2d) then
            kdim = '3D'
            ASSERT(.not. lCell1d)
        end if
!
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=meshNodeCoor)
        if (present(nbCellSupport_)) then
            call utmasu(mesh, kdim, nbSkinCell, cellSkinNume, jvCellVolume, &
                        meshNodeCoor, nbCellSupport_, suppNume_, skinInsideVolume)
        else
            call utmasu(mesh, kdim, nbSkinCell, cellSkinNume, jvCellVolume, &
                        meshNodeCoor, zero, ibid, skinInsideVolume)
        end if
        call jeveuo(jvCellVolume, 'L', vi=vCellVolume)
        do iSkinCell = 1, nbSkinCell
            cellSuppNume(iSkinCell) = vCellVolume(iSkinCell)
        end do
        call jedetr(jvCellVolume)
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkNormalOnSkinCell
!
! Get "volumic" cells support of skin cells
!
! In  mesh             : name of mesh
! In  modelDime        : global dimension of model
! In  nbSkinCell       : number of skin cells
! In  cellSkinNume     : index for skin cells
! In  cellSkinNbNode   : number of nodes for skin cells
! In  cellSkinNodeIndx : index of first node for skin cells
! In  cellSuppNume     : index for "volumic" cells support of skin cells
!
! --------------------------------------------------------------------------------------------------
    subroutine checkNormalOnSkinCell(meshz, modelDime, nbSkinCell, cellSkinNume, cellSkinNbNode, &
                                     cellSkinNodeIndx, cellSuppNume, lMisoriented)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshz
        integer, intent(in) :: modelDime
        integer, intent(in) :: nbSkinCell, cellSkinNume(nbSkinCell), cellSuppNume(nbSkinCell)
        integer, intent(in) :: cellSkinNbNode(nbSkinCell), cellSkinNodeIndx(nbSkinCell)
        aster_logical, intent(out) :: lMisoriented
! - Local
        character(len=8) :: mesh, skinCellType, suppCellType, cellName
        integer :: iSkinCell, suppNume, nodeFirst, iNode
        integer :: skinNume
        integer :: skinNbNode, skinNodeNume(27), suppNbNode, suppNodeNume(27)
        real(kind=8) :: signNorm
        real(kind=8), pointer :: meshNodeCoor(:) => null()
        integer, pointer :: meshConnexLen(:) => null()
        integer, pointer :: meshConnex(:) => null()
        integer, pointer :: meshCellType(:) => null()
!   ------------------------------------------------------------------------------------------------
        mesh = meshz
        call jeveuo(mesh//'.TYPMAIL        ', 'L', vi=meshCellType)
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=meshNodeCoor)
        call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=meshConnexLen)
        call jeveuo(mesh//'.CONNEX', 'E', vi=meshConnex)
        lMisoriented = ASTER_FALSE
!
        do iSkinCell = 1, nbSkinCell
! ----- Support cell of skin
            suppNume = cellSuppNume(iSkinCell)
            if (suppNume .eq. 0) then
                call utmess('A', 'FLUID1_3')
            else
! --------- Parameters of skin cell
                skinNume = cellSkinNume(iSkinCell)
                skinNbNode = cellSkinNbNode(iSkinCell)
                ASSERT(skinNbNode .le. 27)
                nodeFirst = cellSkinNodeIndx(iSKinCell)
                skinNodeNume = 0
                do iNode = 1, skinNbNode
                    skinNodeNume(iNode) = meshConnex(nodeFirst+iNode-1)
                end do
                call jenuno(jexnum('&CATA.TM.NOMTM', meshCellType(skinNume)), skinCellType)
! --------- Parameters of support cell
                nodeFirst = meshConnexLen(suppNume)
                suppNbNode = meshConnexLen(suppNume+1)-nodeFirst
                ASSERT(suppNbNode .le. 27)
                suppNodeNume = 0
                do iNode = 1, suppNbNode
                    suppNodeNume(iNode) = meshConnex(nodeFirst+iNode-1)
                end do
                call jenuno(jexnum('&CATA.TM.NOMTM', meshCellType(suppNume)), suppCellType)
                cellName = int_to_char8(skinNume)
!
! --------- Get sign of normal from skin cell to its volumic support
                call getSignNormalSkinToSupport(modelDime, skinNodeNume, suppNodeNume, &
                                                skinCellType, suppCellType, meshNodeCoor, &
                                                signNorm)
                if (signNorm .gt. 0.d0) then
                    lMisoriented = ASTER_TRUE
                end if
            end if
!
        end do
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getSignNormalSkinToSupport
!
! get sign of normal from skin cell to its volumic support
!
! In  modelDime        : global dimension of model
! In  skinNodeNume     : index of nodes in mesh for skin cell
! In  suppNodeNume     : index of nodes in mesh for support of skin cell
! In  skinCellType     : type of skin cell
! In  suppCellType     : type of support of skin cell
! In  meshNodeCoor     : coordinates of nodes
! Out signNorm         : sign of normal from skin to support
!                          >0 : normal if from skin to inside support cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getSignNormalSkinToSupport(modelDime, skinNodeNume, suppNodeNume, skinCellType, &
                                          suppCellType, meshNodeCoor, signNorm)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer, intent(in) :: modelDime
        integer, intent(in) :: skinNodeNume(27), suppNodeNume(27)
        character(len=8), intent(in) :: skinCellType, suppCellType
        real(kind=8), intent(in) :: meshNodeCoor(*)
        real(kind=8), intent(out) :: signNorm
! - Local
        integer :: node1, node2, node3, iDime, iNode
        integer :: skinNbNode, suppNbNode
        real(kind=8) :: node1Coor(3), node2Coor(3), node3Coor(3)
        real(kind=8) :: n1n3(3), n1n2(3), normal(3), n1g(3), norm
        real(kind=8) :: xgm(3), xg3d(3)
        blas_int :: b_incx, b_incy, b_n
!
!   ------------------------------------------------------------------------------------------------
        signNorm = 0.d0
        normal = 0.d0
!
! - Compute normal to skin element
!
        node1 = skinNodeNume(1)
        node2 = skinNodeNume(2)
        if (modelDime .eq. 3) then
            node3 = skinNodeNume(3)
        end if
        node1Coor = 0.d0
        node2Coor = 0.d0
        do iDime = 1, modelDime
            node1Coor(iDime) = meshNodeCoor(3*(node1-1)+iDime)
            node2Coor(iDime) = meshNodeCoor(3*(node2-1)+iDime)
        end do
        n1n2 = node2Coor-node1Coor
        if (modelDime .eq. 2) then
            normal(1) = n1n2(2)
            normal(2) = -n1n2(1)
        else if (modelDime .eq. 3) then
            do iDime = 1, 3
                node3Coor(iDime) = meshNodeCoor(3*(node3-1)+iDime)
            end do
            n1n3 = node3Coor-node1Coor
            call provec(n1n2, n1n3, normal)
        else
            ASSERT(ASTER_FALSE)
        end if
        call normev(normal, norm)
!
! - Compute barycentric for skin element
!
        if (skinCellType(1:4) .eq. 'QUAD') then
            skinNbNode = 4
        else if (skinCellType(1:4) .eq. 'TRIA') then
            skinNbNode = 3
        else if (skinCellType(1:3) .eq. 'SEG') then
            skinNbNode = 2
            ASSERT(modelDime .eq. 2)
        else
            ASSERT(ASTER_FALSE)
        end if
        xgm = 0.d0
        do iNode = 1, skinNbNode
            xgm(1) = xgm(1)+meshNodeCoor(3*(skinNodeNume(iNode)-1)+1)
            xgm(2) = xgm(2)+meshNodeCoor(3*(skinNodeNume(iNode)-1)+2)
            if (modelDime .eq. 3) then
                xgm(3) = xgm(3)+meshNodeCoor(3*(skinNodeNume(iNode)-1)+3)
            end if
        end do
        xgm(1) = xgm(1)/skinNbNode
        xgm(2) = xgm(2)/skinNbNode
        xgm(3) = xgm(3)/skinNbNode
!
! - Compute barycentric for support element
!
        if (suppCellType(1:4) .eq. 'HEXA') then
            suppNbNode = 8
        else if (suppCellType(1:4) .eq. 'PENT') then
            suppNbNode = 6
        else if (suppCellType(1:4) .eq. 'PYRA') then
            suppNbNode = 5
        else if (suppCellType(1:4) .eq. 'TETR') then
            suppNbNode = 4
        else if (suppCellType(1:4) .eq. 'QUAD') then
            suppNbNode = 4
        else if (suppCellType(1:4) .eq. 'TRIA') then
            suppNbNode = 3
        else
            ASSERT(ASTER_FALSE)
        end if
        xg3d = 0.d0
        do iNode = 1, suppNbNode
            xg3d(1) = xg3d(1)+meshNodeCoor(3*(suppNodeNume(iNode)-1)+1)
            xg3d(2) = xg3d(2)+meshNodeCoor(3*(suppNodeNume(iNode)-1)+2)
            if (modelDime .eq. 3) then
                xg3d(3) = xg3d(3)+meshNodeCoor(3*(suppNodeNume(iNode)-1)+3)
            end if
        end do
        xg3d(1) = xg3d(1)/suppNbNode
        xg3d(2) = xg3d(2)/suppNbNode
        xg3d(3) = xg3d(3)/suppNbNode
!
! - Compute sign of normal from skin to support
!
        n1g = xg3d-xgm
        call normev(n1g, norm)
        b_n = to_blas_int(modelDime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        signNorm = ddot(b_n, n1g, b_incx, normal, b_incy)
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkInclude
!
! check include for mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine checkInclude()
!   ------------------------------------------------------------------------------------------------
! - Local
        integer :: cellTypeNume, nbCell
!   ------------------------------------------------------------------------------------------------
        nbCell = 0
        call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_POI1)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_SEG2)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG3'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_SEG3)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG4'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_SEG4)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA3'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TRIA3)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA6'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TRIA6)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA7'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TRIA7)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD4'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_QUAD4)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD8'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_QUAD8)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD9'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_QUAD9)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TETRA4'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TETRA4)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TETRA10'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TETRA10)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TETRA15'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TETRA15)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA6'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA6)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA15'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA15)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA18'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA18)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA21'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA21)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM5'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PYRAM5)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM13'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PYRAM13)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM19'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PYRAM19)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA8'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_HEXA8)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA20'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_HEXA20)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA27'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_HEXA27)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA9'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_HEXA9)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA7'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA7)
        nbCell = nbCell+1
!
! - Total number of "physical" cells - See MeshTypes_type.h
        ASSERT(nbCell .eq. MT_NPHMAX)
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getPropertiesOfCell
!
! Get properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getPropertiesOfCell(cellTypeName, cellTopo_, cellOrder_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=8), intent(in) :: cellTypeName
        character(len=4), optional, intent(out) :: cellTopo_
        integer, optional, intent(out) :: cellOrder_
! - Local
        integer :: cellOrder
        character(len=4) :: cellTopo
!   ------------------------------------------------------------------------------------------------
        if (cellTypeName .eq. 'POI1') then
            cellTopo = 'PONC'
            cellOrder = 0
        else if (cellTypeName .eq. 'SEG2') then
            cellTopo = 'LINE'
            cellOrder = 1
        else if (cellTypeName .eq. 'SEG3') then
            cellTopo = 'LINE'
            cellOrder = 2
        else if (cellTypeName .eq. 'SEG4') then
            cellTopo = 'LINE'
            cellOrder = 3
        else if (cellTypeName .eq. 'TRIA3') then
            cellTopo = 'SURF'
            cellOrder = 1
        else if (cellTypeName .eq. 'TRIA6') then
            cellTopo = 'SURF'
            cellOrder = 2
        else if (cellTypeName .eq. 'TRIA7') then
            cellTopo = 'SURF'
            cellOrder = 2
        else if (cellTypeName .eq. 'QUAD4') then
            cellTopo = 'SURF'
            cellOrder = 1
        else if (cellTypeName .eq. 'QUAD8') then
            cellTopo = 'SURF'
            cellOrder = 2
        else if (cellTypeName .eq. 'QUAD9') then
            cellTopo = 'SURF'
            cellOrder = 2
        else if (cellTypeName .eq. 'TETRA4') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'TETRA10') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'TETRA15') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PENTA6') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'PENTA15') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PENTA18') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PENTA21') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PYRAM5') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'PYRAM13') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PYRAM19') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'HEXA8') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'HEXA20') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'HEXA27') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'HEXA9') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'PENTA7') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else
            ASSERT(ASTER_FALSE)
        end if
        if (present(cellTopo_)) then
            cellTopo_ = cellTopo
        end if
        if (present(cellOrder_)) then
            cellOrder_ = cellOrder
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCellOptionForName
!
! Get options from user for name of cell
!
! In  keywfact         : factor keyword to read parameters
! In  iocc             : index of factor keyword
! Out lPrefCellName    : flag if a string was given for name of cell
! Out lPrefCellNume    : flag if a integer was given for name of cell
! Out prefCellName     : string for name of cell
! Out prefCellNume     : integer for name of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getCellOptionForName(keywfact, iocc, lPrefCellName, lPrefCellNume, prefCellName, &
                                    prefCellNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: keywfact
        integer, intent(in) :: iocc
        aster_logical, intent(out) :: lPrefCellName, lPrefCellNume
        character(len=8), intent(out) :: prefCellName
        integer, intent(out) :: prefCellNume
! - Local
        integer :: n1
!   ------------------------------------------------------------------------------------------------
        lPrefCellName = ASTER_FALSE
        lPrefCellNume = ASTER_FALSE
        prefCellName = ' '
        prefCellNume = -1
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! createNameOfCell
!
! Create name for cell
!
! IO  cellName         : name of cell
! In  lPrefCellName    : flag if a string was given for name of cell
! In  lPrefCellNume    : flag if a integer was given for name of cell
! In  prefCellName     : string for name of cell
! IO  prefCellNume     : integer for name of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine createNameOfCell(cellName, lPrefCellName, lPrefCellNume, prefCellName, &
                                prefCellNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=8), intent(inout) :: cellName
        aster_logical, intent(in) :: lPrefCellName, lPrefCellNume
        character(len=8), intent(in) :: prefCellName
        integer, intent(inout) :: prefCellNume
! - Local
        integer :: lenPrefName, lenPrefPrev
        character(len=8) :: knume
!   ------------------------------------------------------------------------------------------------
        lenPrefName = lxlgut(prefCellName)
        if (lPrefCellNume) then
            call codent(prefCellNume, 'G', knume)
            prefCellNume = prefCellNume+1
            lenPrefPrev = lxlgut(knume)
            if (lenPrefName+lenPrefPrev .gt. 8) then
                call utmess('F', 'MESH2_1')
            end if
            cellName = knume
        else
            ASSERT(ASTER_FALSE)
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNodeOptionForName
!
! Get options from user for name of node
!
! In  keywfact         : factor keyword to read parameters
! In  iocc             : index of factor keyword
! Out lPrefNodeName    : flag if a string was given for name of node
! Out lPrefNodeNume    : flag if a integer was given for name of node
! Out prefNodeName     : string for name of node
! Out prefNodeNume     : integer for name of node
!
! --------------------------------------------------------------------------------------------------
    subroutine getNodeOptionForName(keywfact, iocc, lPrefNodeName, lPrefNodeNume, prefNodeName, &
                                    prefNodeNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: keywfact
        integer, intent(in) :: iocc
        aster_logical, intent(out) :: lPrefNodeName, lPrefNodeNume
        character(len=8), intent(out) :: prefNodeName
        integer, intent(out) :: prefNodeNume
! - Local
        integer :: n1
!   ------------------------------------------------------------------------------------------------
        lPrefNodeName = ASTER_FALSE
        lPrefNodeNume = ASTER_FALSE
        prefNodeName = ' '
        prefNodeNume = -1
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! createNameOfNode
!
! Create name for node
!
! IO  nodeName         : name of node
! In  lPrefNodeName    : flag if a string was given for name of node
! In  lPrefNodeNume    : flag if a integer was given for name of node
! In  prefNodeName     : string for name of node
! IO  prefNodeNume     : integer for name of node
!
! --------------------------------------------------------------------------------------------------
    subroutine createNameOfNode(nodeName, lPrefNodeName, lPrefNodeNume, prefNodeName, &
                                prefNodeNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=8), intent(inout) :: nodeName
        aster_logical, intent(in) :: lPrefNodeName, lPrefNodeNume
        character(len=8), intent(in) :: prefNodeName
        integer, intent(inout) :: prefNodeNume
! - Local
        integer :: lenPrefName, lenPrefPrev
        character(len=8) :: knume
!   ------------------------------------------------------------------------------------------------
        lenPrefName = lxlgut(prefNodeName)
        if (lPrefNodeNume) then
            call codent(prefNodeNume, 'G', knume)
            prefNodeNume = prefNodeNume+1
            lenPrefPrev = lxlgut(knume)
            if (lenPrefName+lenPrefPrev .gt. 8) then
                call utmess('F', 'MESH2_1')
            end if
            nodeName = knume
        else
            ASSERT(ASTER_FALSE)
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getGroupsFromCell
!
! Get some groups where cells are
!
! In  mesh             : mesh
! In  cellNume         : index of current cell
! Out groupCell        : name of groups of cells
! Out nbGroupCell      : number of group
!
! --------------------------------------------------------------------------------------------------
    subroutine getGroupsFromCell(meshZ, cellNume, groupCell, nbGroupCell)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        integer, intent(in) :: cellNume
        character(len=24), intent(out) :: groupCell(4)
        integer, intent(out) :: nbGroupCell
! - Local
        character(len=8) :: mesh
        integer :: iexi, nbGroup, iGroup, iCell, nbCell
        character(len=24) :: groupName
        integer, pointer :: cellList(:) => null()
!   ------------------------------------------------------------------------------------------------
        mesh = meshZ
        nbGroupCell = 0
        groupCell = " "
        call jeexin(mesh//'.GROUPEMA', iexi)
        if (iexi .gt. 0) then
            call jelira(mesh//'.GROUPEMA', 'NMAXOC', nbGroup)
            do iGroup = 1, nbGroup
                call jeexin(jexnum(mesh//'.GROUPEMA', iGroup), iexi)
                if (iexi .eq. 0) cycle
                call jenuno(jexnum(mesh//'.GROUPEMA', iGroup), groupName)
                call jelira(jexnum(mesh//'.GROUPEMA', iGroup), 'LONUTI', nbCell)
                call jeveuo(jexnum(mesh//'.GROUPEMA', iGroup), 'L', vi=cellList)
                do iCell = 1, nbCell
                    if (cellList(iCell) .eq. cellNume) then
                        nbGroupCell = nbGroupCell+1
                        groupCell(nbGroupCell) = groupName
                        if (nbGroupCell .eq. 4) goto 100
                    end if
                end do
            end do
        end if
100     continue
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMeshDimension
!
! Get dimension of mesh
!
! In  mesh             : mesh
! Out meshDime         : dimension of mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine getMeshDimension(meshZ, meshDime)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ
        integer, intent(out) :: meshDime
! ----- Local
        character(len=8) :: mesh
        character(len=24) :: answer
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        meshDime = 0
        call dismoi('Z_CST', mesh, 'MAILLAGE', repk=answer)
        if (answer .eq. 'OUI') then
            meshDime = 2
        else
            meshDime = 3
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getFirstNodeFromNodeGroup
!
! Get first node from group of nodes
!
! In  mesh             : mesh
! In  nodeGroup        ; name of group of node
! Out nodeName         : name of first node from group of nodes
! Out nodeNume         : index of first node from group of nodes
!
! --------------------------------------------------------------------------------------------------
    subroutine getFirstNodeFromNodeGroup(meshZ, nodeGroupZ, nodeName, nodeNume)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ, nodeGroupZ
        character(len=8), intent(out) :: nodeName
        integer, intent(out) :: nodeNume
! ----- Local
        character(len=8) :: mesh
        character(len=24) :: nodeGroup
        integer :: iret
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        nodeGroup = nodeGroupZ
        nodeName = " "
        nodeNume = 0
        call utnono(" ", mesh, 'NOEUD', nodeGroup, nodeName, &
                    iret)
        if (iret .eq. 0) then
            nodeNume = char8_to_int(nodeName)
        else if (iret .eq. 10) then
            call utmess('F', 'MESH3_1', sk=nodeGroup)
        else if (iret .eq. 1) then
            call utmess('A', 'MESH3_2', sk=nodeGroup)
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListOfCellGroup
!
! Get list of groups of cells
!
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine getListOfCellGroup(meshZ, factorKeywordZ, iFactorKeyword, nbGroup, listOfGroups)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ, factorKeywordZ
        integer, intent(in) :: iFactorKeyword
        character(len=24), pointer :: listOfGroups(:)
        integer, intent(out) :: nbGroup
! ----- Local
        character(len=8) :: mesh
        character(len=24) :: k24Dummy
        integer :: nbRet
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        nbGroup = 0
        call getvem(mesh, 'GROUP_MA', factorKeywordZ, 'GROUP_MA', iFactorKeyword, &
                    0, k24Dummy, nbGroup)
        nbGroup = -nbGroup
        if (nbGroup .ne. 0) then
            allocate (listOfGroups(nbGroup))
            call getvem(mesh, 'GROUP_MA', factorKeywordZ, 'GROUP_MA', iFactorKeyword, &
                        nbGroup, listOfGroups, nbRet)
            ASSERT(nbRet .eq. nbGroup)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkCellsAreSkin
!
! Check if cells are skins' ones
!
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine checkCellsAreSkin(meshZ, nbCell, listCellNume, onlySkin1D, listCellType, &
                                 hasSkin1D, hasSkin2D)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ
        integer, intent(in) :: nbCell
        integer, pointer :: listCellNume(:)
        aster_logical, intent(in) :: onlySkin1D
        character(len=8), pointer :: listCellType(:)
        aster_logical, intent(out) :: hasSkin1D, hasSkin2D
! ----- Local
        character(len=8) :: mesh
        integer :: cellNume, iCell, cellTypeNume
        character(len=8) :: cellTypeName
        integer, pointer :: typmail(:) => null()
!
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        hasSkin1D = ASTER_FALSE
        hasSkin2D = ASTER_FALSE
!
! ----- Access to mesh datastructures
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
!
! ----- Check cells
        do iCell = 1, nbCell
            cellNume = listCellNume(iCell)
!
! --------- Get type of cell
            cellTypeNume = typmail(cellNume)
            call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
            listCellType(iCell) = cellTypeName
!
! --------- Detect type
            if (cellTypeName(1:4) .eq. 'QUAD') then
                hasSkin2D = ASTER_TRUE
            else if (cellTypeName(1:4) .eq. 'TRIA') then
                hasSkin2D = ASTER_TRUE
            else if (cellTypeName(1:3) .eq. 'SEG') then
                hasSkin1D = ASTER_TRUE
            else
                call utmess('F', 'MESH3_94', sk=cellTypeName)
            end if
            if (hasSkin1D .and. hasSkin2D) then
                call utmess('F', 'MESH3_98')
            end if
        end do
!
        if (onlySkin1D .and. hasSkin2D) then
            call utmess('F', 'MESH3_92')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module mesh_module
