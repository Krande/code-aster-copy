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
module SolidShell_Mesh_module
! ==================================================================================================
use crea_maillage_module
! ==================================================================================================
implicit none
! ==================================================================================================
public  :: orieHexa9, setValueOnFace
private :: isThisQuad
! ==================================================================================================
private 
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getelem.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/jexnum.h"
#include "asterfort/jenuno.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! orieHexa9
!
! Orientation of Hexa9 cells for solid shell
!
! In  iOcc             : index of factor keyword COQUE_SOLIDE
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
subroutine orieHexa9(iOcc, nbVoluCell, voluCell, mesh)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    integer, intent(in) :: iOcc, nbVoluCell
    integer, pointer :: voluCell(:)
    character(len=8), intent(in) :: mesh
! - Local
    character(len=16), parameter :: keywfact  = 'COQUE_SOLIDE'
    integer :: iVoluNode, iVoluCell, iSurfCell, iSurfNode
    integer :: voluCellNume, surfCellNume
    integer :: nbSurfCell, voluNbNode, surfNbNode, voluCellType, surfCellType
    integer :: voluNodeHexa(8), surfNodeQuad(4)
    character(len=16) :: suffix, answer
    character(len=24), parameter :: jvSurfCell = '&&ORIHEXA9.SURF'
    integer :: fp(24), iter, aux
    integer, pointer :: typmail(:) => null()
    integer, pointer :: surfCell(:) => null()
    integer, pointer :: surfNode(:) => null()
    integer, pointer :: voluNode(:) => null()
    aster_logical :: lFace1, lFace2, lHexaInMesh, lPentaInMesh
!   ------------------------------------------------------------------------------------------------
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
    call dismoi('EXI_HEXA8', mesh, 'MAILLAGE', repk = answer)
    lHexaInMesh  = answer .eq. 'OUI'
    call dismoi('EXI_PENTA6', mesh, 'MAILLAGE', repk = answer)
    lPentaInMesh = answer .eq. 'OUI'

! - Get list of surfacic cells
    suffix    = '_SURF'
    call getelem(mesh      , keywfact, iocc, ' ', jvSurfCell,&
                 nbSurfCell, suffix)
    if (nbSurfCell .eq. 0) then
        if (lHexaInMesh) then
            call utmess('F', 'SOLIDSHELL1_2')
        endif
        goto 99
    endif
    if (lPentaInMesh .and. (.not. lHexaInMesh)) then
        call utmess('A', 'SOLIDSHELL1_3')
    endif
    call jeveuo(jvSurfCell, 'L', vi = surfCell)

! - Reorient
    do iVoluCell = 1, nbVoluCell
        voluCellNume = voluCell(iVoluCell)
        voluCellType = typmail(voluCellNume)
        if (voluCellType .ne. MT_HEXA8) then
            call utmess('F', 'SOLIDSHELL1_1')
        endif
        call jelira(jexnum(mesh//'.CONNEX', voluCellNume), 'LONMAX', voluNbNode)
        call jeveuo(jexnum(mesh//'.CONNEX', voluCellNume), 'E', vi = voluNode)
        do iSurfCell = 1, nbSurfCell
            surfCellNume = surfCell(iSurfCell)
            surfCellType = typmail(surfCellNume)
            if (surfCellType .ne. MT_QUAD4) then
                call utmess('F', 'SOLIDSHELL1_1')
            endif
            call jelira(jexnum(mesh//'.CONNEX', surfCellNume), 'LONMAX', surfNbNode)
            call jeveuo(jexnum(mesh//'.CONNEX', surfCellNume), 'L', vi = surfNode)
            do iVoluNode = 1, 8
                voluNodeHexa(iVoluNode) = voluNode(iVoluNode)
            end do
            
            do iSurfNode = 1, 4
                surfNodeQuad(iSurfNode) = surfNode(iSurfNode)
            end do
            fp(1)  = voluNodeHexa(1)
            fp(2)  = voluNodeHexa(2)
            fp(3)  = voluNodeHexa(3)
            fp(4)  = voluNodeHexa(4)
            fp(5)  = voluNodeHexa(5)
            fp(6)  = voluNodeHexa(6)
            fp(7)  = voluNodeHexa(7)
            fp(8)  = voluNodeHexa(8)
            fp(9)  = voluNodeHexa(5)
            fp(10) = voluNodeHexa(1)
            fp(11) = voluNodeHexa(4)
            fp(12) = voluNodeHexa(8)
            fp(13) = voluNodeHexa(6)
            fp(14) = voluNodeHexa(2)
            fp(15) = voluNodeHexa(3)
            fp(16) = voluNodeHexa(7)
            fp(17) = voluNodeHexa(5)
            fp(18) = voluNodeHexa(1)
            fp(19) = voluNodeHexa(2)
            fp(20) = voluNodeHexa(6)
            fp(21) = voluNodeHexa(8)
            fp(22) = voluNodeHexa(4)
            fp(23) = voluNodeHexa(3)
            fp(24) = voluNodeHexa(7)
            do iter = 1, 3
                aux = 8*(iter-1)
                lFace1 = isThisQuad(fp(1+aux:4+aux), surfNodeQuad(:))
                lFace2 = isThisQuad(fp(5+aux:8+aux), surfNodeQuad(:))
                if (lFace1 .or. lFace2) then
                    voluNode(1) = fp(1+aux)
                    voluNode(2) = fp(2+aux)
                    voluNode(3) = fp(3+aux)
                    voluNode(4) = fp(4+aux)
                    voluNode(5) = fp(5+aux)
                    voluNode(6) = fp(6+aux)
                    voluNode(7) = fp(7+aux)
                    voluNode(8) = fp(8+aux)
                end if
            end do
        end do
    end do
!
99  continue
!
    call jedetr(jvSurfCell)
!
!   ------------------------------------------------------------------------------------------------
end subroutine
! --------------------------------------------------------------------------------------------------
!
! isThisQuad
!
! Check if quad is the good one
!
! In  quadToTest       : list of nodes of quad to test
! In  quadRefe         : list of nodes of quad
!
! --------------------------------------------------------------------------------------------------
function isThisQuad(quadToTest, quadRefe)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    aster_logical :: isThisQuad
    integer, intent(in), dimension(:) :: quadToTest, quadRefe
! - Local
    integer :: ii, jj, nbNodeFound
! --------------------------------------------------------------------------------------------------
!
    isThisQuad = ASTER_FALSE
    nbNodeFound  = 0
    do ii = 1, 4
        do jj = 1, 4
            if (quadToTest(ii) .eq. quadRefe(jj)) then
                nbNodeFound = nbNodeFound + 1
            endif
        end do
    end do
    isThisQuad = nbNodeFound .eq. 4
!
!   ------------------------------------------------------------------------------------------------
end function
! --------------------------------------------------------------------------------------------------
!
! setValueOnFace
!
! Set value on face of volumic element
! From list of skin element, find volumic cell underlying and set value (real or string) on
! inferior face or superior face of this volumic cell
!
! In  mesh             : mesh
! In  valeR            : value to set on face (real)
! In  valeK            : value to set on face (string)
! In  nbCellVolu       : total number of volumic cells
! Ptr cellVolu         : pointer to volumic cells
! In  nbCellSkin       : total number of skin cells
! Ptr cellSkin         : pointer to skin cells
! Ptr valeCell         : pointer to index of volumic cell where value is affected
! Ptr valeFaceR        : pointer to value (real) on face
! Ptr valeFaceK        : pointer to value (string) on face
!
! --------------------------------------------------------------------------------------------------
subroutine setValueOnFace(mesh      ,&
                          valeR     , valeK   ,&
                          nbCellVolu, cellVolu,&
                          nbCellSkin, cellSkin,&
                          valeCell  , valeFaceR , valeFaceK)
!   ------------------------------------------------------------------------------------------------
! - Parameters
    character(len=8), intent(in) :: mesh
    real(kind=8)    , intent(in) :: valeR
    character(len=8), intent(in) :: valeK
    integer, intent(in)          :: nbCellVolu
    integer, pointer             :: cellVolu(:)
    integer, intent(in)          :: nbCellSkin
    integer, pointer             :: cellSkin(:)
    integer, pointer             :: valeCell(:)
    real(kind=8), pointer        :: valeFaceR(:)
    character(len=8), pointer    :: valeFaceK(:)
! - Local
    integer :: iCellSkin, iNode, iCellVolu
    integer :: cellSkinNume, cellVoluNume, cellTypeNume
    character(len=8) :: cellTypeName, cellSkinName
    integer :: nodeSkin(4), nodeFaceInf(4), nodeFaceSup(4)
    integer, pointer :: meshConnex(:) => null(), meshTypmail(:) => null()
    aster_logical :: lFaceInf, lFaceSup
!   ------------------------------------------------------------------------------------------------
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi = meshTypmail)
!
    do iCellSkin = 1, nbCellSkin

! ----- Get current face
        cellSkinNume = cellSkin(iCellSkin)
        ASSERT(cellSkinNume .gt. 0)
        cellTypeNume = meshTypmail(cellSkinNume)
        call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
        if (cellTypeName .ne. 'QUAD4') then
            call jenuno(jexnum(mesh//'.NOMMAI', cellSkinNume), cellSkinName)
            call utmess('F', 'CHARGES_9', sk = cellSkinName)
        endif

! ----- Get nodes
        call jeveuo(jexnum(mesh//'.CONNEX', cellSkinNume), 'E', vi = meshConnex)
        do iNode = 1, 4
            nodeSkin(iNode) = meshConnex(iNode)
        end do

! ----- Search volume and set to faces of this volume
        lFaceInf = ASTER_FALSE
        lFaceSup = ASTER_FALSE
        do iCellVolu = 1, nbCellVolu
            cellVoluNume = cellVolu(iCellVolu)
            call jeveuo(jexnum(mesh//'.CONNEX', cellVoluNume), 'E', vi = meshConnex)
            do iNode = 1, 4
                nodeFaceInf(iNode) = meshConnex(iNode)
                nodeFaceSup(iNode) = meshConnex(iNode+4)
            end do
            lFaceInf = isThisQuad(nodeSkin, nodeFaceInf)
            valeCell(iCellVolu) = cellVoluNume
            if (lFaceInf) then
                valeFaceR(2*(iCellVolu-1)+1) = valeR
                valeFaceK(2*(iCellVolu-1)+1) = valeK
            endif
            lFaceSup = isThisQuad(nodeSkin, nodeFaceSup)
            if (lFaceSup) then
                valeFaceR(2*(iCellVolu-1)+2) = valeR
                valeFaceK(2*(iCellVolu-1)+2) = valeK
            endif
            if (lFaceInf .or. lFaceSup) then
                exit
            endif
        end do
        if (.not. lFaceInf .and. .not. lFaceSup) then
            call utmess('F', 'CHARGES_10')
        endif
        if (lFaceInf .and. lFaceSup) then
            ASSERT(ASTER_FALSE)
        endif
    end do
!
!   ------------------------------------------------------------------------------------------------
end subroutine
!
end module SolidShell_Mesh_module
