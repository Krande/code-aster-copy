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
!
! ==================================================================================================
!
! Module for elementary computation for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Elementary_module
! ==================================================================================================
use SolidShell_type
use SolidShell_Utilities_module
use SolidShell_Debug_module
use SolidShell_Geometry_module
use SolidShell_Elementary_Hexa_module
! ==================================================================================================
implicit none
! ==================================================================================================
public  :: compRigiMatr
private :: setMateOrientation, compElemElasMatrix,&
           initCellGeom, initMateProp, initElemProp
! ==================================================================================================
private
#include "jeveux.h"
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/SolidShell_type.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/dmat3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcangm.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! initElemProp
!
! Initialization of general properties of finite element
!
! In  inteFami         : name of integration scheme family
! Out elemProp         : general properties of element
!
! --------------------------------------------------------------------------------------------------
subroutine initElemProp(inteFami, elemProp)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    character(len=4), intent(in)     :: inteFami
    type(SSH_ELEM_PROP), intent(out) :: elemProp
! - Local
    integer :: npg, nno
    integer :: jvWeight, jvCoor, jvShape, jvDShape
!   ------------------------------------------------------------------------------------------------
!
    if (SSH_DBG_ELEM) SSH_DBG_STRG('> initElemProp')

! - Get element parameters
    call elrefe_info(fami = inteFami, npg = npg, nno = nno,&
                     jpoids = jvWeight, jcoopg = jvCoor,&
                     jvf    = jvShape , jdfde=jvDShape)

! - Set parameters for integration scheme
    elemProp%elemInte%inteFami    = inteFami
    elemProp%elemInte%nbIntePoint = npg
    elemProp%elemInte%jvCoor      = jvCoor
    elemProp%elemInte%jvWeight    = jvWeight
    elemProp%elemInte%jvShape     = jvShape
    elemProp%elemInte%jvDShape    = jvDShape

! - Set main parameters of finite element
    if (nno .eq. 9) then
        elemProp%cellType = SSH_CELL_HEXA
    else
        ASSERT(ASTER_FALSE)
    endif
    elemProp%nbNode     = nno
    elemProp%nbNodeGeom = nno - 1
    elemProp%nbDofGeom  = 3*elemProp%nbNodeGeom
    elemProp%nbDof      = elemProp%nbDofGeom + 1
!
    if (SSH_DBG_ELEM) SSH_DBG_STRG('< initElemProp')
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! initCellGeom
!
! Initialization of geometric properties of cell
!
! In  elemProp         : general properties of element
! Out cellGeom         : general geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
subroutine initCellGeom(elemProp, cellGeom)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in)  :: elemProp
    type(SSH_CELL_GEOM), intent(out) :: cellGeom
! - Local
    integer      :: iDofGeom, iNodeGeom
    real(kind=8) :: detJac0
!   ------------------------------------------------------------------------------------------------
!
    if (SSH_DBG_ELEM) SSH_DBG_STRG('> initCellGeom')

! - Access to field of coordinates
    call jevech('PGEOMER', 'L', cellGeom%jvGeom)

! - Set initial geometry
    do iDofGeom = 1, elemProp%nbDofGeom
        cellGeom%geomInit(iDofGeom) = zr(cellGeom%jvGeom+iDofGeom-1)
    end do
    do iNodeGeom = 1, elemProp%nbNodeGeom
        cellGeom%geomInitX(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+1)
        cellGeom%geomInitY(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+2)
        cellGeom%geomInitZ(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+3)
    enddo

! - Set center
    if (elemProp%cellType == SSH_CELL_HEXA) then
        cellGeom%cellCenterCova = hexaCovaCenter
    else
        ASSERT(ASTER_FALSE)
    endif

! - Compute Jacobian matrix at center of element on initial configuration
    call compJacoMatr(elemProp         ,&
                      cellGeom%geomInit, cellGeom%cellCenterCova,&
                      cellGeom%Jac0    , cellGeom%JacInv0       ,&
                      detJac0)
    cellGeom%detJac0 = abs(detJac0)
!
    if (SSH_DBG_ELEM) SSH_DBG_STRG('< initCellGeom')
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! initMateProp
!
! Initialization of properties of material
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! Out matePara         : parameters of material
!
! --------------------------------------------------------------------------------------------------
subroutine initMateProp(elemProp, cellGeom, matePara)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in)  :: elemProp
    type(SSH_CELL_GEOM), intent(in)  :: cellGeom
    type(SSH_MATE_PARA), intent(out) :: matePara
! - Local
    integer :: jvMate
!   ------------------------------------------------------------------------------------------------
!
    if (SSH_DBG_ELEM) SSH_DBG_STRG('> initMateProp')
! - Access to field of material parameters
    call jevech('PMATERC', 'L', jvMate)
    matePara%jvMater = zi(jvMate)

! - Set material orientation
    call setMateOrientation(elemProp, cellGeom, matePara)

! - Compute elasticity matrix at middle of cell
    call compElemElasMatrix(elemProp%elemInte, matePara)

!
    if (SSH_DBG_ELEM) SSH_DBG_STRG('< initMateProp')
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! setMateOrientation
!
! Set material orientation
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  matePara         : parameters of material
!
! --------------------------------------------------------------------------------------------------
subroutine setMateOrientation(elemProp, cellGeom, matePara)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_PROP), intent(in)    :: elemProp
    type(SSH_CELL_GEOM), intent(in)    :: cellGeom
    type(SSH_MATE_PARA), intent(inout) :: matePara
! - Local
    real(kind=8) :: bary(3)
    integer :: iNode, iDime, nno, jvGeom
!   ------------------------------------------------------------------------------------------------
!
    nno    = elemProp%nbNodeGeom
    jvGeom = cellGeom%jvGeom

! - Compute barycentric center
    bary = 0.d0
    do iNode = 1, nno
        do iDime = 1, SSH_NDIM
            bary(iDime) = bary(iDime) + zr(jvGeom+iDime+SSH_NDIM*(iNode-1)-1)/nno
        end do
    end do

! - Get orientation
    call rcangm(SSH_NDIM, bary, matePara%mateBase)
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compElemElasMatrix
!
! Compute elasticity matrix at middle of cell
!
! In  elemInte         : properties of integration scheme
! IO  matePara         : parameters of material
!
! --------------------------------------------------------------------------------------------------
subroutine compElemElasMatrix(elemInte, matePara)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    type(SSH_ELEM_INTE), intent(in)    :: elemInte
    type(SSH_MATE_PARA), intent(inout) :: matePara
! - Local
    real(kind=8) :: xyzgau(3)
!   ------------------------------------------------------------------------------------------------
! 
   xyzgau = 0.d0
   call dmat3d(elemInte%inteFami, matePara%jvMater , r8vide(), '+', 1,&
               1                , matePara%mateBase, xyzgau  , matePara%elemHookeMatrix)
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! compRigiMatr
!
! Compute rigidity matrix - RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
subroutine compRigiMatr()
!   ------------------------------------------------------------------------------------------------
! - Local
    character(len=4), parameter :: inteFami = 'RIGI'
    type(SSH_CELL_GEOM) :: cellGeom
    type(SSH_ELEM_PROP) :: elemProp
    type(SSH_MATE_PARA) :: matePara
    real(kind=8) :: matrRigi(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
    integer :: jvMatr, i, j, k
!   ------------------------------------------------------------------------------------------------
!
    matrRigi = 0.d0

! - Initialization of general properties of finite element
    call initElemProp(inteFami, elemProp)
    if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! - Initialization of geometric properties of cell
    call initCellGeom(elemProp, cellGeom)
    if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! - Initialization of properties of material
    call initMateProp(elemProp, cellGeom, matePara)
    if (SSH_DBG_MATE) call dbgObjMatePara(matePara)

! - Compute rigidity matrix
    if (elemProp%cellType .eq. SSH_CELL_HEXA) then
        call compRigiMatrHexa(elemProp, cellGeom, matePara, matrRigi)
    else
        ASSERT(ASTER_FALSE)
    endif

! - Save matrix
    call jevech('PMATUUR', 'E', jvMatr)
    k = 0
    do i = 1, elemProp%nbDof
        do j = 1, i
            k = k + 1
            zr(jvMatr-1+k) = matrRigi(i, j)
        end do
    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
!
end module SolidShell_Elementary_module
