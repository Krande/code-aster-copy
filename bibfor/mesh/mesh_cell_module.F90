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
! ==================================================================================================
!
! Module for cells in mesh
!
! ==================================================================================================
!
module mesh_cell_module
! ==================================================================================================
    use mesh_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: cellCreate
    public :: cellCompTang, cellCompNorm, cellCompBaseAtPoint, cellCompNormAtBary
    public :: cellCompDiam, cellCompBary, cellCompCenterGlob
    public :: cellDebug, cellCopyType, cellPoinParaToGlob
    public :: cellPoinInside, cellPoinAdjust
    public :: cellSetType
    private :: cellSetCoorGlob, cellSetLine, cellSetUndef
    private :: cellCompCenterPara, cellSetNbNeigh
    private :: cellPoinAdjustSeg, cellPoinAdjustTria, cellPoinAdjustQuad
! ==================================================================================================
    private
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfvf.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! cellCreate
!
! Create a cell
!
! In  cellNume         : index of cell in mesh
! Ptr nodeCoor         : pointer to coordinates of nodes
! Ptr meshTypeGeom     : pointer to type of cells in mesh
! Ptr meshConx         : pointers to connectivity of mesh
!     meshConxCumu
! Out cellGeom         : general geometric properties of cell
! Out cellGeomLine     : general geometric properties of cell after linearization
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCreate(cellNume, nodeCoor, &
                          meshTypeGeom, meshConx, meshConxCumu, &
                          cellGeom, cellGeomLine_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: cellNume
        real(kind=8), pointer :: nodeCoor(:)
        integer(kind=8), pointer :: meshTypeGeom(:), meshConx(:), meshConxCumu(:)
        type(CELL_GEOM), intent(out) :: cellGeom
        type(CELL_GEOM), optional, intent(out) :: cellGeomLine_
!   ------------------------------------------------------------------------------------------------
!

! ----- Set type of cell from mesh
        call cellSetType(meshTypeGeom, cellNume, cellGeom)

! ----- Set coordinates of cell (in global space)
        call cellSetCoorGlob(cellNume, nodeCoor, &
                             meshConx, meshConxCumu, &
                             cellGeom)

! ----- Compute barycenters
        call cellCompBary(cellGeom)

! ----- Compute diameter
        call cellCompDiam(cellGeom)

! ----- Linearize cell
        if (present(cellGeomLine_)) then
            call cellSetLine(cellGeom, cellGeomLine_)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellSetType
!
! Get type of cell
!
! Ptr meshTypeGeom     : geometric type of cells
! In  cellNume         : index of cell in mesh
! IO  cellGeom         : geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellSetType(meshTypeGeom, cellNume, cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), pointer :: meshTypeGeom(:)
        integer(kind=8), intent(in) :: cellNume
        type(CELL_GEOM), intent(inout) :: cellGeom
! ----- Local
        integer(kind=8) :: cellTypeNume, nbNode, nbNodeS, cellDime
        character(len=8) :: cellTypeName
!   ------------------------------------------------------------------------------------------------
!
        cellTypeNume = meshTypeGeom(cellNume)
        call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
        select case (cellTypeName)
        case ('SEG2')
            cellGeom%cellCode = 'SE2'
        case ('SEG3')
            cellGeom%cellCode = 'SE3'
        case ('TRIA3')
            cellGeom%cellCode = 'TR3'
        case ('TRIA6')
            cellGeom%cellCode = 'TR6'
        case ('QUAD4')
            cellGeom%cellCode = 'QU4'
        case ('QUAD8')
            cellGeom%cellCode = 'QU8'
        case ('QUAD9')
            cellGeom%cellCode = 'QU9'
        case default
            ASSERT(ASTER_FALSE)
        end select

! ----- Prov: all cells are skin's ones
        cellGeom%isSkin = ASTER_TRUE

! ----- For standard cells
        call elrfno(cellGeom%cellCode, nbNode, nbNodeS, cellDime, nodeCoor=cellGeom%coorNodePara)
        cellGeom%nbNode = nbNode
        cellGeom%cellDime = cellDime
        cellGeom%isLinear = (nbNode .eq. nbNodeS)
        if (cellGeom%isSkin) then
            call cellSetNbNeigh(cellGeom)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellSetCoorGlob
!
! Set coordinates of cell (in global reference frame)
!
! In  cellNume         : index of cell in mesh
! Ptr nodeCoor         : pointer to coordinates of nodes for this cell
! Ptr meshConnex       : pointer to connectivity of mesh
! Ptr meshConnexCumu   : pointer to connectivity of mesh
! IO  cellGeom         : geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellSetCoorGlob(cellNume, nodeCoor, &
                               meshConnex, meshConnexCumu, &
                               cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: cellNume
        real(kind=8), pointer :: nodeCoor(:)
        integer(kind=8), pointer :: meshConnex(:), meshConnexCumu(:)
        type(CELL_GEOM), intent(inout) :: cellGeom
! ----- Local
        integer(kind=8) :: iNode, iDime, nodeNume
!   ------------------------------------------------------------------------------------------------
!
        cellGeom%coorNodeGlob = 0.d0
        do iNode = 1, cellGeom%nbNode
            nodeNume = meshConnex(meshConnexCumu(cellNume)-1+iNode)
            do iDime = 1, 3
                cellGeom%coorNodeGlob(iDime, iNode) = &
                    nodeCoor(3*(nodeNume-1)+iDime)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellSetLine
!
! Linearization of cell
!
! In  cellGeom         : geometric properties of cell
! Out cellGeomLine     : geometry of cell after linearization
!
! --------------------------------------------------------------------------------------------------
    subroutine cellSetLine(cellGeom, cellGeomLine)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        type(CELL_GEOM), intent(out) :: cellGeomLine
!   ------------------------------------------------------------------------------------------------
!
        cellGeomLine = cellGeom

! ----- Change support
        if (cellGeom%cellCode .eq. 'SE2' .or. cellGeom%cellCode .eq. 'SE3') then
            cellGeomLine%cellCode = 'SE2'
            cellGeomLine%nbNode = 2
        elseif (cellGeom%cellCode .eq. 'TR3' .or. cellGeom%cellCode .eq. 'TR6') then
            cellGeomLine%cellCode = 'TR3'
            cellGeomLine%nbNode = 3
        else if (cellGeom%cellCode .eq. 'QU4' .or. &
                 cellGeom%cellCode .eq. 'QU8' .or. &
                 cellGeom%cellCode .eq. 'QU9') then
            cellGeomLine%cellCode = 'QU4'
            cellGeomLine%nbNode = 4
        else
            WRITE (6, *) "cellGeom%cellCode: ", cellGeom%cellCode
            ASSERT(ASTER_FALSE)
        end if
        cellGeomLine%isLinear = ASTER_TRUE

! ----- Change coordinates
        cellGeomLine%coorNodeGlob(:, cellGeomLine%nbNode+1:MT_NNOMAX3D) = 0.d0
        cellGeomLine%coorNodePara(:, cellGeomLine%nbNode+1:MT_NNOMAX3D) = 0.d0

! ----- Update barycenter
        call cellCompBary(cellGeomLine)

! ----- Update diameter
        call cellCompDiam(cellGeomLine)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellSetUndef
!
! Set undefined values for cell
!
! IO  cellGeom         : geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellSetUndef(cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(inout) :: cellGeom
!   ------------------------------------------------------------------------------------------------
!
        cellGeom%diameter = r8nnem()
        cellGeom%baryGlob = r8nnem()
        cellGeom%baryPara = r8nnem()
        cellGeom%coorNodeGlob = r8nnem()
        cellGeom%coorNodePara = r8nnem()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCompTang
!
! Compute tangents on cell at current point
!
! In  cellGeom         : geometric properties of cell
! In  spaceDime        : dimension of space (2 or 3)
! In  dShapeFunc       : values of derivatives of shape functions at current point
! Out tau1             : first tangent at current point
! Out tau2             : second tangent at current point
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCompTang(cellGeom, spaceDime, dShapeFunc, tau1, tau2)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        integer(kind=8), intent(in) :: spaceDime
        real(kind=8), intent(in) :: dShapeFunc(3, MT_NNOMAX3D)
        real(kind=8), intent(out) :: tau1(3), tau2(3)
! ----- Local
        real(kind=8), parameter :: zero = 0.d0
        integer(kind=8) :: iNode, iDime
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(cellGeom%nbNode .le. MT_NNOMAX3D)
        tau1 = zero
        tau2 = zero
        do iDime = 1, spaceDime
            do iNode = 1, cellGeom%nbNode
                tau1(iDime) = cellGeom%coorNodeGlob(iDime, iNode)*dShapeFunc(1, iNode)+tau1(iDime)
                if (spaceDime .eq. 3) then
                    tau2(iDime) = cellGeom%coorNodeGlob(iDime, iNode)*dShapeFunc(2, iNode)+ &
                                  tau2(iDime)
                end if
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCompBaseAtPoint
!
! Compute base at point (outward normal)
!
! In  cellGeom         : geometric properties of cell
! In  spaceDime        : dimension of space (2 or 3)
! In  poinCoorPara     : coordinates of point in parametric space
! Out baseExte         : base on point (with outward normal)
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCompBaseAtPoint(cellGeom, spaceDime, poinCoorPara, baseExte)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        integer(kind=8), intent(in) :: spaceDime
        real(kind=8), intent(in) :: poinCoorPara(3)
        type(CELL_SKIN_BASE), intent(out) :: baseExte
! ----- Local
        real(kind=8) :: dShapeFunc(3, MT_NNOMAX3D), tau1(3), tau2(3), norm(3)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(cellGeom%isSkin)
        baseExte%spaceDime = spaceDime
        baseExte%normIsExte = ASTER_TRUE
        baseExte%norm = 0.d0
        baseExte%tau = 0.d0
        call elrfdf(cellGeom%cellCode, poinCoorPara, dShapeFunc)
        call cellCompTang(cellGeom, spaceDime, dShapeFunc, tau1, tau2)
        call cellCompNorm(spaceDime, tau1, tau2, norm)
        baseExte%norm = -norm
        baseExte%tau(:, 1) = tau1
        baseExte%tau(:, 2) = tau2
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCompNorm
!
! Compute norm
!
! In  spaceDime        : dimension of space (2 or 3)
! In  tau1             : first tangent
! In  tau2             : second tangent
! Out norm             : norm
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCompNorm(spaceDime, tau1, tau2, norm)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: spaceDime
        real(kind=8), intent(in) :: tau1(3), tau2(3)
        real(kind=8), intent(out) :: norm(3)
! ----- Local
        real(kind=8) :: noor
!   ------------------------------------------------------------------------------------------------
!
        norm = 0.d0
        noor = 0.d0
        if (spaceDime .eq. 2) then
            norm(1) = -tau1(2)
            norm(2) = tau1(1)
            norm(3) = 0.d0
        else if (spaceDime .eq. 3) then
            call provec(tau2, tau1, norm)
        else
            ASSERT(ASTER_FALSE)
        end if
        call normev(norm, noor)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellPoinParaToGlob
!
! Change coordinates of a point in cell to global reference frame
!
! In  cellGeom         : geometric properties of cell
! In  poinCoorPara     : coordinates of point in cell parametric frame
! Out poinCoorGlob     : coordinates of point in global reference frame
!
! --------------------------------------------------------------------------------------------------
    subroutine cellPoinParaToGlob(cellGeom, &
                                  poinCoorPara, poinCoorGlob)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in) :: poinCoorPara(2)
        real(kind=8), intent(out) :: poinCoorGlob(3)
! ----- Local
        real(kind=8) :: shapeFunc(MT_NNOMAX3D)
        integer(kind=8) :: iNode, iDime
!   ------------------------------------------------------------------------------------------------
!
        poinCoorGlob = 0.d0

! ----- Get shape functions in parametric space
        shapeFunc = 0.d0
        call elrfvf(cellGeom%cellCode, poinCoorPara, shapeFunc)

! ----- Change coordinates
        do iDime = 1, 3
            do iNode = 1, cellGeom%nbNode
                poinCoorGlob(iDime) = poinCoorGlob(iDime)+ &
                                      cellGeom%coorNodeGlob(iDime, iNode)*shapeFunc(iNode)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCompCenterPara
!
! Compute center of cell in parametric space
!
! In  cellGeom         : geometric properties of cell
! Out cellCentPara     : center of cell in parametric space
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCompCenterPara(cellGeom, cellCentPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(out) :: cellCentPara(3)
! ----- Local
        integer(kind=8) :: iNode
!   ------------------------------------------------------------------------------------------------
!
        cellCentPara = 0.d0
        do iNode = 1, cellGeom%nbNode
            cellCentPara = cellCentPara+cellGeom%coorNodePara(:, iNode)
        end do
        cellCentPara = cellCentPara/cellGeom%nbNode
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCompCenterGlob
!
! Compute center of cell in global space
!
! In  cellGeom         : geometric properties of cell
! Out cellCentPara     : center of cell in parametric space
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCompCenterGlob(cellGeom, cellCentGlob)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(out) :: cellCentGlob(3)
! ----- Local
        real(kind=8) :: cellCentPara(3)
!   ------------------------------------------------------------------------------------------------
!
        cellCentGlob = 0.d0
        call cellCompCenterPara(cellGeom, cellCentPara)
        call cellPoinParaToGlob(cellGeom, cellCentPara, cellCentGlob)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCompDiam
!
! Compute diameter of cell
!
! IO  cellGeom         : geometry of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCompDiam(cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(inout) :: cellGeom
! ----- Locals
        integer(kind=8) :: iNode, jNode
        real(kind=8) :: length
!   ------------------------------------------------------------------------------------------------
!
        cellGeom%diameter = 0.d0
        do iNode = 1, cellGeom%nbNode
            do jNode = iNode+1, cellGeom%nbNode
                length = norm2(cellGeom%coorNodeGlob(1:3, iNode)- &
                               cellGeom%coorNodeGlob(1:3, jNode))
                cellGeom%diameter = max(cellGeom%diameter, length)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCompBary
!
! Compute barycenters of cell
!
! In  spaceDime        : dimension of space (2 or 3)
! IO  cellGeom         : geometry of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCompBary(cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(inout) :: cellGeom
! ----- Local
        real(kind=8) :: cellBaryPara(3)
        integer(kind=8) :: iNode
!   ------------------------------------------------------------------------------------------------
!
        cellGeom%baryGlob = 0.d0
        do iNode = 1, cellGeom%nbNode
            cellGeom%baryGlob = cellGeom%baryGlob+ &
                                cellGeom%coorNodeGlob(:, iNode)
        end do
        cellGeom%baryGlob = cellGeom%baryGlob/real(cellGeom%nbNode, kind=8)
        call cellCompCenterPara(cellGeom, cellBaryPara)
        cellGeom%baryPara = cellBaryPara
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCompNormAtBary
!
! Compute normal at barycenter of cell
!
! In  spaceDime        : dimension of space (2 or 3)
! In  cellGeom         : geometric properties of cell
! Out cellNormBary     : normal at barycenter of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCompNormAtBary(spaceDime, cellGeom, cellNormBary)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: spaceDime
        type(CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(out) :: cellNormBary(3)
! ----- Local
        real(kind=8) :: cellCentPara(3)
        type(CELL_SKIN_BASE) :: baseExte
!   ------------------------------------------------------------------------------------------------
!
        cellNormBary = 0.d0
        ASSERT(cellGeom%isSkin)

! ----- Get center of cell in parametric space
        cellCentPara = cellGeom%baryPara

! ----- Compute base on point (outward normal)
        call cellCompBaseAtPoint(cellGeom, spaceDime, cellCentPara, baseExte)
        cellNormBary = baseExte%norm
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellSetNbNeigh
!
! Set number of neighbours of skin cell
!
! IO  cellGeom         : geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellSetNbNeigh(cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(inout) :: cellGeom
! ----- Local
        integer(kind=8) :: nbNeigh
!   ------------------------------------------------------------------------------------------------
!
        nbNeigh = 0
        ASSERT(cellGeom%isSkin)
        if (cellGeom%cellCode == 'SE2' .or. &
            cellGeom%cellCode == 'SE3') then
            nbNeigh = 2
        elseif (cellGeom%cellCode == 'TR3' .or. &
                cellGeom%cellCode == 'TR6' .or. &
                cellGeom%cellCode == 'TR7') then
            nbNeigh = 3
        elseif (cellGeom%cellCode == 'QU4' .or. &
                cellGeom%cellCode == 'QU8' .or. &
                cellGeom%cellCode == 'QU9') then
            nbNeigh = 4
        else
            ASSERT(ASTER_FALSE)
        end if
        cellGeom%nbNeigh = nbNeigh
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellDebug
!
! Print debug informations about cell
!
! In  cellGeom         : geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellDebug(cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
! ----- Locals
        integer(kind=8) :: iNode
!   ------------------------------------------------------------------------------------------------
!
        if (cellGeom%isLinear) then
            WRITE (6, *) "Linéaire"
        else
            WRITE (6, *) "Quadratique"
        end if
        if (cellGeom%isSkin) then
            WRITE (6, *) "Maille de bord"
        else
            WRITE (6, *) "Maille d'intérieur"
        end if
        WRITE (6, *) "Type      :", cellGeom%cellCode
        WRITE (6, *) "Dime      :", cellGeom%cellDime
        WRITE (6, *) "Nb node   :", cellGeom%nbNode
        do iNode = 1, cellGeom%nbNode
            WRITE (6, *) "Node (", iNode, "): ", &
                cellGeom%coorNodeGlob(1:3, iNode), " (global frame)"
            WRITE (6, *) "Node (", iNode, "): ", &
                cellGeom%coorNodePara(1:3, iNode), " (local frame)"
        end do
        WRITE (6, *) "Diameter  :", cellGeom%diameter
        WRITE (6, *) "Barycenter:", cellGeom%baryGlob, " (global frame)"
        WRITE (6, *) "Barycenter:", cellGeom%baryPara, " (local frame)"
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellCopyType
!
! Create new cell with same geometric type
!
! In  cellGeom         : geometric properties of original cell
! Out cellCopy         : geometric properties of new cell
!
! --------------------------------------------------------------------------------------------------
    subroutine cellCopyType(cellOrig, cellCopy)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellOrig
        type(CELL_GEOM), intent(out) :: cellCopy
! ----- Locals

!   ------------------------------------------------------------------------------------------------
!
        cellCopy = cellOrig
        call cellSetUndef(cellCopy)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine

! --------------------------------------------------------------------------------------------------
!
! cellPoinInside
!
! Detect if point is inside cell - With its own reference frame
!
! In  toleInside       : tolerance to detect
! In  cellGeom         : geometric properties of cell
! In  poinCoorPara     : parametric coordinates of point
! Out poinIsInside
!
! --------------------------------------------------------------------------------------------------
    function cellPoinInside(cellGeom, toleInside, poinCoorPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: cellPoinInside
        type(CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in) :: toleInside
        real(kind=8), intent(in) :: poinCoorPara(2)
! ----- Local
        character(len=8) :: cellCode
!   ------------------------------------------------------------------------------------------------
!
        cellCode = cellGeom%cellCode
        ASSERT(cellGeom%isSkin)
        cellPoinInside = ASTER_FALSE
        if (cellCode .eq. 'SE2' .or. cellCode .eq. 'SE3') then
            if (poinCoorPara(1) .ge. (-1.d0-toleInside) .and. &
                poinCoorPara(1) .le. (1.d0+toleInside)) then
                cellPoinInside = ASTER_TRUE
            end if
        elseif (cellCode .eq. 'TR3' .or. cellCode .eq. 'TR6') then
            if (poinCoorPara(1) .ge. -toleInside .and. &
                poinCoorPara(2) .ge. -toleInside .and. &
                (poinCoorPara(2)+poinCoorPara(1)) .le. (1.d0+toleInside)) then
                cellPoinInside = ASTER_TRUE
            end if
        elseif (cellCode .eq. 'QU4' .or. cellCode .eq. 'QU8' .or. cellCode .eq. 'QU9') then
            if (poinCoorPara(1) .ge. (-1.d0-toleInside) .and. &
                poinCoorPara(1) .le. (1.d0+toleInside) .and. &
                poinCoorPara(2) .ge. (-1.d0-toleInside) .and. &
                poinCoorPara(2) .le. (1.d0+toleInside)) then
                cellPoinInside = ASTER_TRUE
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! cellPoinAdjust
!
! Adjust parametric coordinates to be inside cell
!
! In  cellGeom         : geometric properties of cell
! In  projOutside      : tolerance for projection outside cell
! IO  poinCoorPara     : coordinates of point
! Out projType         : type of projection
!                        0 - No projection inside cell (point is inside)
!                        1 - Projection inside cell
!                        2 - No projection inside cell (to far away)
!
! --------------------------------------------------------------------------------------------------
    subroutine cellPoinAdjust(cellGeom, projOutside, &
                              poinCoorPara, projType_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in) :: projOutside
        real(kind=8), intent(inout) :: poinCoorPara(2)
        integer(kind=8), optional, intent(out) :: projType_
! ----- Local
        integer(kind=8) :: projType
        character(len=8) :: cellCode
!   ------------------------------------------------------------------------------------------------
!
        projType = 0
        cellCode = cellGeom%cellCode
        ASSERT(cellGeom%isSkin)

        if (cellCode(1:2) .eq. 'SE') then
            call cellPoinAdjustSeg(projOutside, &
                                   poinCoorPara, projType)
        else if (cellCode(1:2) .eq. 'TR') then
            call cellPoinAdjustTria(projOutside, &
                                    poinCoorPara, projType)
        else if (cellCode(1:2) .eq. 'QU') then
            call cellPoinAdjustQuad(projOutside, &
                                    poinCoorPara, projType)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (present(projType_)) then
            projType_ = projType
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellPoinAdjustSeg
!
! Adjust parametric coordinates to be inside cell when cell is segment
!
! In  projOutside      : tolerance for projection outside cell
! IO  poinCoorPara     : coordinates of point
! Out projType         : type of projection
!                        0 - No projection inside cell (point is inside)
!                        1 - Projection inside cell
!                        2 - No projection inside cell (to far away)
!
! --------------------------------------------------------------------------------------------------
    subroutine cellPoinAdjustSeg(projOutside, &
                                 poinCoorPara, projType)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: projOutside
        real(kind=8), intent(inout) :: poinCoorPara(2)
        integer(kind=8), intent(out) :: projType
! ----- Local
        aster_logical :: near
        real(kind=8) :: ecart
!   tolerances --- absolue et relative --- pour determiner si deux distances sont egales
        real(kind=8), parameter :: atol = 1.e-12
        real(kind=8), parameter :: rtol = 1.e-12
!   ------------------------------------------------------------------------------------------------
!
        projType = 0
        ecart = -1.d0

! ----- Premier ajustement : on positionne le point sur le bord, s'il est à une distance
! ----- (normalisée) inférieure a atol du bord
        if (abs(poinCoorPara(1)+1.d0) .le. atol) then
            poinCoorPara(1) = -1.d0
        end if
        if (abs(poinCoorPara(1)-1.d0) .le. atol) then
            poinCoorPara(1) = +1.d0
        end if

! ----- RABATTEMENT
        if ((poinCoorPara(1) .lt. -1.d0) .or. (poinCoorPara(1) .gt. 1.d0)) then
            ecart = abs(poinCoorPara(1))-1.d0
            projType = 1
            if (poinCoorPara(1) .lt. -1.d0) then
                poinCoorPara(1) = -1.d0
            else if (poinCoorPara(1) .gt. 1.d0) then
                poinCoorPara(1) = 1.d0
            end if
            near = abs(ecart-projOutside) .le. (atol+projOutside*rtol)
            if (ecart .gt. projOutside .and. .not. near) then
                projType = 2
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellPoinAdjustTria
!
! Adjust parametric coordinates to be inside cell when cell is triangle
!
! In  projOutside      : tolerance for projection outside cell
! IO  poinCoorPara     : coordinates of point
! Out projType         : type of projection
!                        0 - No projection inside cell (point is inside)
!                        1 - Projection inside cell
!                        2 - No projection inside cell (to far away)
!
! --------------------------------------------------------------------------------------------------
    subroutine cellPoinAdjustTria(projOutside, &
                                  poinCoorPara, projType)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: projOutside
        real(kind=8), intent(inout) :: poinCoorPara(2)
        integer(kind=8), intent(out) :: projType
! ----- Local
        integer(kind=8) :: izone
        aster_logical :: near
        real(kind=8) :: ecart, k1pk2, k2mk1, ksi1e, ksi2e
!   tolerances --- absolue et relative --- pour determiner si deux distances sont egales
        real(kind=8), parameter :: atol = 1.e-12
        real(kind=8), parameter :: rtol = 1.e-12
!   ------------------------------------------------------------------------------------------------
!
        projType = 0
        ecart = -1.d0
        k1pk2 = poinCoorPara(1)+poinCoorPara(2)
        k2mk1 = poinCoorPara(2)-poinCoorPara(1)

! ----- Premier ajustement : on positionne le point sur le bord, s'il est à une distance
! ----- (normalisée) inférieure a atol du bord
        if (abs(poinCoorPara(1)) .le. atol) then
            poinCoorPara(1) = 0.d0
        end if
        if (abs(poinCoorPara(2)) .le. atol) then
            poinCoorPara(2) = 0.d0
        end if
        if (abs(k1pk2-1.d0) .le. atol) then
            k1pk2 = +1.d0
        end if
        if (abs(k2mk1+1.d0) .le. atol) then
            k2mk1 = -1.d0
        end if
        if (abs(k2mk1-1.d0) .le. atol) then
            k2mk1 = +1.d0
        end if
        if ((poinCoorPara(1) .ge. 0.d0) .and. (poinCoorPara(2) .ge. 0.d0) .and. &
            (k1pk2 .le. 1.d0)) then
            goto 99
        end if

! ----- SECTEUR CONCERNE
        izone = 0
        if (poinCoorPara(1) .lt. 0.d0) then
            if (poinCoorPara(2) .lt. 0.d0) then
                izone = 1
            else if ((poinCoorPara(2) .ge. 0.d0) .and. (poinCoorPara(2) .le. 1.d0)) then
                izone = 2
            else if (poinCoorPara(2) .gt. 1.d0) then
                izone = 3
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if (poinCoorPara(2) .lt. 0.d0) then
            if (poinCoorPara(1) .lt. 0.d0) then
                izone = 1
            else if ((poinCoorPara(1) .ge. 0.d0) .and. (poinCoorPara(1) .le. 1.d0)) then
                izone = 8
            else if (poinCoorPara(1) .gt. 1.d0) then
                izone = 7
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if (poinCoorPara(1) .ge. 0.d0) then
            if (k2mk1 .gt. 1.d0) then
                izone = 4
            elseif ((k1pk2 .gt. 1.d0) .and. (k2mk1 .ge. -1.d0) &
                    .and. (k2mk1 .le. 1.d0)) then
                izone = 5
                ksi1e = 5.d-1*(1.d0+poinCoorPara(1)-poinCoorPara(2))
                ksi2e = 5.d-1*(1.d0-poinCoorPara(1)+poinCoorPara(2))
            else if ((poinCoorPara(2) .ge. 0.d0) .and. (k2mk1 .lt. -1.d0)) then
                izone = 6
            end if
        end if

! ----- CALCUL DE L'ECART
        if (izone .eq. 1) then
            ecart = sqrt(abs(poinCoorPara(1))*abs(poinCoorPara(1))+ &
                         abs(poinCoorPara(2))*abs(poinCoorPara(2)))
        else if (izone .eq. 2) then
            ecart = sqrt(abs(poinCoorPara(1))*abs(poinCoorPara(1)))
        else if (izone .eq. 3 .or. izone .eq. 4) then
            ecart = sqrt(abs(poinCoorPara(1))*abs(poinCoorPara(1))+ &
                         (poinCoorPara(2)-1.d0)*(poinCoorPara(2)-1.d0))
        else if (izone .eq. 5) then
            ecart = sqrt((poinCoorPara(1)-ksi1e)*(poinCoorPara(1)-ksi1e)+ &
                         (poinCoorPara(2)-ksi2e)*(poinCoorPara(2)-ksi2e))
        else if (izone .eq. 6 .or. izone .eq. 7) then
            ecart = sqrt(abs(poinCoorPara(2))*abs(poinCoorPara(2))+ &
                         (poinCoorPara(1)-1.d0)*(poinCoorPara(1)-1.d0))
        else if (izone .eq. 8) then
            ecart = sqrt(abs(poinCoorPara(2))*abs(poinCoorPara(2)))
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- RABATTEMENT
        projType = 1
        if (izone .eq. 1) then
            poinCoorPara(1) = 0.d0
            poinCoorPara(2) = 0.d0
        else if (izone .eq. 2) then
            poinCoorPara(1) = 0.d0
        else if (izone .eq. 3 .or. izone .eq. 4) then
            poinCoorPara(1) = 0.d0
            poinCoorPara(2) = 1.d0
        else if (izone .eq. 5) then
            poinCoorPara(1) = ksi1e
            poinCoorPara(2) = ksi2e
        else if (izone .eq. 6 .or. izone .eq. 7) then
            poinCoorPara(1) = 1.d0
            poinCoorPara(2) = 0.d0
        else if (izone .eq. 8) then
            poinCoorPara(2) = 0.d0
        end if

        near = abs(ecart-projOutside) .le. (atol+projOutside*rtol)
        if (ecart .gt. projOutside .and. .not. near) then
            projType = 2
        end if
99      continue
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! cellPoinAdjustQuad
!
! Adjust parametric coordinates to be inside cell when cell is quadrangle
!
! In  projOutside      : tolerance for projection outside cell
! IO  poinCoorPara     : coordinates of point
! Out projType         : type of projection
!                        0 - No projection inside cell (point is inside)
!                        1 - Projection inside cell
!                        2 - No projection inside cell (to far away)
!
! --------------------------------------------------------------------------------------------------
    subroutine cellPoinAdjustQuad(projOutside, &
                                  poinCoorPara, projType)
!   ------------------------------------------------------------------------------------------------
! ----- Parameter
        real(kind=8), intent(in) :: projOutside
        real(kind=8), intent(inout) :: poinCoorPara(2)
        integer(kind=8), intent(out) :: projType
! ----- Local
        integer(kind=8) :: izone
        aster_logical :: near
        real(kind=8) :: ecart, k1pk2, k2mk1
!   tolerances --- absolue et relative --- pour determiner si deux distances sont egales
        real(kind=8), parameter :: atol = 1.e-12
        real(kind=8), parameter :: rtol = 1.e-12
!   ------------------------------------------------------------------------------------------------
!
        projType = 0
        ecart = -1.d0
        k1pk2 = poinCoorPara(1)+poinCoorPara(2)
        k2mk1 = poinCoorPara(2)-poinCoorPara(1)

! ----- Premier ajustement : on positionne le point sur le bord, s'il est à une distance
! ----- (normalisée) inférieure a atol du bord
        if (abs(poinCoorPara(1)+1.d0) .le. atol) then
            poinCoorPara(1) = -1.d0
        end if
        if (abs(poinCoorPara(1)-1.d0) .le. atol) then
            poinCoorPara(1) = +1.d0
        end if
        if (abs(poinCoorPara(2)+1.d0) .le. atol) then
            poinCoorPara(2) = -1.d0
        end if
        if (abs(poinCoorPara(2)-1.d0) .le. atol) then
            poinCoorPara(2) = +1.d0
        end if
        if ((poinCoorPara(1) .ge. -1.d0) .and. (poinCoorPara(1) .le. 1.d0) .and. &
            (poinCoorPara(2) .ge. -1.d0) .and. (poinCoorPara(2) .le. 1.d0)) then
            goto 99
        end if

! ----- SECTEUR CONCERNE
        izone = 0
        if (poinCoorPara(1) .lt. -1.d0) then
            if (poinCoorPara(2) .lt. -1.d0) then
                izone = 1
            else if ((poinCoorPara(2) .ge. -1.d0) .and. (poinCoorPara(2) .le. 1.d0)) then
                izone = 2
            else if (poinCoorPara(2) .gt. 1.d0) then
                izone = 3
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if (poinCoorPara(1) .gt. 1.d0) then
            if (poinCoorPara(2) .lt. -1.d0) then
                izone = 7
            else if ((poinCoorPara(2) .ge. -1.d0) .and. (poinCoorPara(2) .le. 1.d0)) then
                izone = 6
            else if (poinCoorPara(2) .gt. 1.d0) then
                izone = 5
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
        if ((poinCoorPara(1) .ge. -1.d0) .and. (poinCoorPara(1) .le. 1.d0)) then
            if (poinCoorPara(2) .lt. -1.d0) then
                izone = 8
            else if (poinCoorPara(2) .gt. 1.d0) then
                izone = 4
            end if
        end if

! ----- CALCUL DE L'ECART
        if (izone .eq. 1) then
            ecart = sqrt((abs(poinCoorPara(1))-1.d0)*(abs(poinCoorPara(1))-1.d0)+ &
                         (abs(poinCoorPara(2))-1.d0)*(abs(poinCoorPara(2))-1.d0))
        else if (izone .eq. 2) then
            ecart = sqrt((abs(poinCoorPara(1))-1.d0)*(abs(poinCoorPara(1))-1.d0))
        else if (izone .eq. 3) then
            ecart = sqrt((abs(poinCoorPara(1))-1.d0)*(abs(poinCoorPara(1))-1.d0)+ &
                         (abs(poinCoorPara(2))-1.d0)*(abs(poinCoorPara(2))-1.d0))
        else if (izone .eq. 4) then
            ecart = sqrt((abs(poinCoorPara(2))-1.d0)*(abs(poinCoorPara(2))-1.d0))
        else if (izone .eq. 5) then
            ecart = sqrt((abs(poinCoorPara(1))-1.d0)*(abs(poinCoorPara(1))-1.d0)+ &
                         (abs(poinCoorPara(2))-1.d0)*(abs(poinCoorPara(2))-1.d0))
        else if (izone .eq. 6) then
            ecart = sqrt((abs(poinCoorPara(1))-1.d0)*(abs(poinCoorPara(1))-1.d0))
        else if (izone .eq. 7) then
            ecart = sqrt((abs(poinCoorPara(1))-1.d0)*(abs(poinCoorPara(1))-1.d0)+ &
                         (abs(poinCoorPara(2))-1.d0)*(abs(poinCoorPara(2))-1.d0))
        else if (izone .eq. 8) then
            ecart = sqrt((abs(poinCoorPara(2))-1.d0)*(abs(poinCoorPara(2))-1.d0))
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- RABATTEMENT
        projType = 1
        if (izone .eq. 1) then
            poinCoorPara(1) = -1.d0
            poinCoorPara(2) = -1.d0
        else if (izone .eq. 2) then
            poinCoorPara(1) = -1.d0
        else if (izone .eq. 3) then
            poinCoorPara(1) = -1.d0
            poinCoorPara(2) = 1.d0
        else if (izone .eq. 4) then
            poinCoorPara(2) = 1.d0
        else if (izone .eq. 5) then
            poinCoorPara(1) = 1.d0
            poinCoorPara(2) = 1.d0
        else if (izone .eq. 6) then
            poinCoorPara(1) = 1.d0
        else if (izone .eq. 7) then
            poinCoorPara(1) = 1.d0
            poinCoorPara(2) = -1.d0
        else if (izone .eq. 8) then
            poinCoorPara(2) = -1.d0
        end if

        near = abs(ecart-projOutside) .le. (atol+projOutside*rtol)
        if (ecart .gt. projOutside .and. .not. near) then
            projType = 2
        end if
!
99      continue
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module
