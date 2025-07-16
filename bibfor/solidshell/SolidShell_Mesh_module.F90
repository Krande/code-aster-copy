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
module SolidShell_Mesh_module
! ==================================================================================================
    use crea_maillage_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public  :: orieHexa9, setValueOnFace
    private :: isThisQuad, getVolumeLinkedToSurface
! ==================================================================================================
    private
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/cncinv.h"
#include "asterfort/dismoi.h"
#include "asterfort/getelem.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/jeexin.h"
#include "asterfort/int_to_char8.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getVolumeLinkedToSurface
!
! Get volumic cell (give by voluCellType) linked to surface cell (always QUAD4)
!
! --------------------------------------------------------------------------------------------------
    subroutine getVolumeLinkedToSurface(typmail, cnxinvLoncum, cnxinvCell, &
                                        cellSkinNode, voluCellType, &
                                        cellVoluNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), pointer :: typmail(:)
        integer(kind=8), pointer :: cnxinvLoncum(:), cnxinvCell(:)
        integer(kind=8), pointer :: cellSkinNode(:)
        integer(kind=8), intent(in) :: voluCellType
        integer(kind=8), intent(out) :: cellVoluNume
! - Local
        integer(kind=8) :: nbVoluLinkedToSurf, iVoluLinked
        integer(kind=8) :: nbCellLinkedToNode, iCellLinked, cellLinkedNume
        integer(kind=8) :: iNodeSkin, nodeSkinCurr
        integer(kind=8) :: cellLinkedType
        integer(kind=8), parameter :: nbMaxConnected = 15
        integer(kind=8) :: voluLinkedToNode(nbMaxConnected, 2)
        aster_logical :: lVoluFind
!   ------------------------------------------------------------------------------------------------
!
        nbVoluLinkedToSurf = 0
        voluLinkedToNode = 0
        cellVoluNume = 0

        do iNodeSkin = 1, 4
            nodeSkinCurr = cellSkinNode(iNodeSkin)
            nbCellLinkedToNode = &
                cnxinvLoncum(nodeSkinCurr+1)-cnxinvLoncum(nodeSkinCurr)
            do iCellLinked = 1, nbCellLinkedToNode
                cellLinkedNume = cnxinvCell(cnxinvLoncum(nodeSkinCurr)-1+iCellLinked)
                cellLinkedType = typmail(cellLinkedNume)
                if (cellLinkedType .eq. voluCellType) then
                    lVoluFind = ASTER_FALSE
                    do iVoluLinked = 1, nbVoluLinkedToSurf
                        if (voluLinkedToNode(iVoluLinked, 1) .eq. cellLinkedNume) then
                            voluLinkedToNode(iVoluLinked, 2) = &
                                voluLinkedToNode(iVoluLinked, 2)+1
                            lVoluFind = ASTER_TRUE
                        end if
                    end do
                    if (.not. lVoluFInd) then
                        nbVoluLinkedToSurf = nbVoluLinkedToSurf+1
                        ASSERT(nbVoluLinkedToSurf .le. nbMaxConnected)
                        voluLinkedToNode(nbVoluLinkedToSurf, 1) = cellLinkedNume
                        voluLinkedToNode(nbVoluLinkedToSurf, 2) = 1
                    end if
                end if
            end do
        end do

! - Find the volume cell with all nodes belong to surface cell
        do iVoluLinked = 1, nbVoluLinkedToSurf
            if (voluLinkedToNode(iVoluLinked, 2) .eq. 4) then
                cellVoluNume = voluLinkedToNode(iVoluLinked, 1)
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
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
    subroutine orieHexa9(iOcc, mesh)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: iOcc
        character(len=8), intent(in) :: mesh
! - Local
        character(len=16), parameter :: keywfact = 'COQUE_SOLIDE'
        integer(kind=8) :: iCellSkin
        integer(kind=8) :: cellVoluNume, cellSkinNume
        integer(kind=8) :: nbCellSkin, cellSkinType
        integer(kind=8) :: iExist, iDUmmy
        character(len=16) :: suffix, answer
        character(len=19), parameter :: cnxinv = '&&ORIHEXA9.INV'
        character(len=24), parameter :: jvCellSkin = '&&ORIHEXA9.SURF'
        integer(kind=8) :: fp(24), iter, aux
        integer(kind=8), pointer :: typmail(:) => null()
        integer(kind=8), pointer :: surfCell(:) => null()
        integer(kind=8), pointer :: cellVoluNode(:) => null(), cellSkinNode(:) => null()
        integer(kind=8), pointer :: cnxinvLoncum(:) => null(), cnxinvCell(:) => null()
        aster_logical :: lCellSkin1, lCellSkin2, lHexaInMesh, lPentaInMesh
!   ------------------------------------------------------------------------------------------------
!
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
        call dismoi('EXI_HEXA8', mesh, 'MAILLAGE', repk=answer)
        lHexaInMesh = answer .eq. 'OUI'
        call dismoi('EXI_PENTA6', mesh, 'MAILLAGE', repk=answer)
        lPentaInMesh = answer .eq. 'OUI'

! - Prepare reverse connectivity
        call jeexin(cnxinv, iExist)
        if (iExist .eq. 0) then
            iDummy = 0
            call cncinv(mesh, [iDummy], 0, 'V', cnxinv)
        end if
        call jeveuo(jexatr(cnxinv, 'LONCUM'), 'L', vi=cnxinvLoncum)
        call jeveuo(jexnum(cnxinv, 1), 'L', vi=cnxinvCell)

! - Get list of skin cells
        suffix = '_SURF'
        call getelem(mesh, keywfact, iocc, ' ', jvCellSkin, &
                     nbCellSkin, suffix)
        if (nbCellSkin .eq. 0) then
            if (lHexaInMesh) then
                call utmess('F', 'SOLIDSHELL1_2')
            end if
            goto 99
        end if
        if (lPentaInMesh .and. (.not. lHexaInMesh)) then
            call utmess('A', 'SOLIDSHELL1_3')
        end if

! - Loop on skin cells
        call jeveuo(jvCellSkin, 'L', vi=surfCell)
        do iCellSkin = 1, nbCellSkin
            cellSkinNume = surfCell(iCellSkin)
            cellSkinType = typmail(cellSkinNume)
            if (cellSkinType .ne. MT_QUAD4) then
                call utmess('F', 'SOLIDSHELL1_1')
            end if
            ASSERT(cellSkinNume .ne. 0)

! ----- Look for volumic cells attached to current surfacic cell
            call jeveuo(jexnum(mesh//'.CONNEX', cellSkinNume), 'L', vi=cellSkinNode)
            call getVolumeLinkedToSurface(typmail, cnxinvLoncum, cnxinvCell, &
                                          cellSkinNode, MT_HEXA8, &
                                          cellVoluNume)
            ASSERT(cellVoluNume .ne. 0)

! ----- Access to volumic cell
            call jeveuo(jexnum(mesh//'.CONNEX', cellVoluNume), 'E', vi=cellVoluNode)
            fp(1) = cellVoluNode(1)
            fp(2) = cellVoluNode(2)
            fp(3) = cellVoluNode(3)
            fp(4) = cellVoluNode(4)
            fp(5) = cellVoluNode(5)
            fp(6) = cellVoluNode(6)
            fp(7) = cellVoluNode(7)
            fp(8) = cellVoluNode(8)
            fp(9) = cellVoluNode(5)
            fp(10) = cellVoluNode(1)
            fp(11) = cellVoluNode(4)
            fp(12) = cellVoluNode(8)
            fp(13) = cellVoluNode(6)
            fp(14) = cellVoluNode(2)
            fp(15) = cellVoluNode(3)
            fp(16) = cellVoluNode(7)
            fp(17) = cellVoluNode(5)
            fp(18) = cellVoluNode(1)
            fp(19) = cellVoluNode(2)
            fp(20) = cellVoluNode(6)
            fp(21) = cellVoluNode(8)
            fp(22) = cellVoluNode(4)
            fp(23) = cellVoluNode(3)
            fp(24) = cellVoluNode(7)
            do iter = 1, 3
                aux = 8*(iter-1)
                lCellSkin1 = isThisQuad(fp(1+aux:4+aux), cellSkinNode(:))
                lCellSkin2 = isThisQuad(fp(5+aux:8+aux), cellSkinNode(:))
                if (lCellSkin1 .or. lCellSkin2) then
                    cellVoluNode(1) = fp(1+aux)
                    cellVoluNode(2) = fp(2+aux)
                    cellVoluNode(3) = fp(3+aux)
                    cellVoluNode(4) = fp(4+aux)
                    cellVoluNode(5) = fp(5+aux)
                    cellVoluNode(6) = fp(6+aux)
                    cellVoluNode(7) = fp(7+aux)
                    cellVoluNode(8) = fp(8+aux)
                end if
            end do
        end do
!
99      continue
!
        call jedetr(jvCellSkin)
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
        integer(kind=8), intent(in), dimension(:) :: quadToTest, quadRefe
! - Local
        integer(kind=8) :: ii, jj, nbNodeFound
! --------------------------------------------------------------------------------------------------
!
        isThisQuad = ASTER_FALSE
        nbNodeFound = 0
        do ii = 1, 4
            do jj = 1, 4
                if (quadToTest(ii) .eq. quadRefe(jj)) then
                    nbNodeFound = nbNodeFound+1
                end if
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
! In  nbCellSkin       : total number of skin cells
! Ptr cellSkin         : pointer to skin cells
! Ptr valeCell         : pointer to index of volumic cell where value is affected
! Ptr valeFaceR        : pointer to value (real) on face
! Ptr valeFaceK        : pointer to value (string) on face
!
! --------------------------------------------------------------------------------------------------
    subroutine setValueOnFace(mesh, &
                              valeR, valeK, &
                              nbCellSkin, cellSkin, &
                              valeCell, valeFaceR, valeFaceK)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=8), intent(in) :: mesh
        real(kind=8), intent(in) :: valeR
        character(len=8), intent(in) :: valeK
        integer(kind=8), intent(in)          :: nbCellSkin
        integer(kind=8), pointer             :: cellSkin(:)
        integer(kind=8), pointer             :: valeCell(:)
        real(kind=8), pointer        :: valeFaceR(:)
        character(len=8), pointer    :: valeFaceK(:)
! - Local
        integer(kind=8) :: iCellSkin, iCellVolu
        integer(kind=8) :: cellSkinNume, cellVoluNume, cellSkinType
        character(len=8) :: cellSkinName
        integer(kind=8) :: voluSkinNodeInf(4), voluSkinNodeSup(4)
        integer(kind=8), pointer :: cellSkinNode(:) => null(), meshTypmail(:) => null()
        integer(kind=8), pointer :: cellVoluNode(:) => null()
        aster_logical :: lVoluSkinInf, lVoluSkinSup
        character(len=19), parameter :: cnxinv = '&&ORIHEXA9.INV'
        integer(kind=8) :: iExist, iDummy

        integer(kind=8), pointer :: cnxinvLoncum(:) => null(), cnxinvCell(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=meshTypmail)

! - Prepare reverse connectivity
        call jeexin(cnxinv, iExist)
        if (iExist .eq. 0) then
            iDummy = 0
            call cncinv(mesh, [iDummy], 0, 'V', cnxinv)
        end if
        call jeveuo(jexatr(cnxinv, 'LONCUM'), 'L', vi=cnxinvLoncum)
        call jeveuo(jexnum(cnxinv, 1), 'L', vi=cnxinvCell)
!
        iCellVolu = 0
        do iCellSkin = 1, nbCellSkin

! ----- Get current skin cell
            cellSkinNume = cellSkin(iCellSkin)
            ASSERT(cellSkinNume .gt. 0)
            cellSkinType = meshTypmail(cellSkinNume)
            if (cellSkinType .ne. MT_QUAD4) then
                cellSkinName = int_to_char8(cellSkinNume)
                call utmess('F', 'SOLIDSHELL1_7', sk=cellSkinName)
            end if

! ----- Get nodes of this skin cell
            call jeveuo(jexnum(mesh//'.CONNEX', cellSkinNume), 'E', vi=cellSkinNode)

! ----- Look for volumic cells attached to current surfacic cell
            call getVolumeLinkedToSurface(meshTypmail, cnxinvLoncum, cnxinvCell, &
                                          cellSkinNode, MT_HEXA9, &
                                          cellVoluNume)
            if (cellVoluNume .eq. 0) then
                call utmess('F', 'SOLIDSHELL1_6')
            end if

! ----- Get nodes of this volume
            call jeveuo(jexnum(mesh//'.CONNEX', cellVoluNume), 'E', vi=cellVoluNode)

! ----- Search superior or inferior face
            voluSkinNodeInf(1:4) = cellVoluNode(1:4)
            voluSkinNodeSup(1:4) = cellVoluNode(5:8)
            lVoluSkinInf = isThisQuad(cellSkinNode, voluSkinNodeInf)
            valeCell(iCellSkin) = cellVoluNume
            if (lVoluSkinInf) then
                valeFaceR(2*(iCellSkin-1)+1) = valeR
                valeFaceK(2*(iCellSkin-1)+1) = valeK
            end if
            lVoluSkinSup = isThisQuad(cellSkinNode, voluSkinNodeSup)
            if (lVoluSkinSup) then
                valeFaceR(2*(iCellSkin-1)+2) = valeR
                valeFaceK(2*(iCellSkin-1)+2) = valeK
            end if
            if (.not. lVoluSkinInf .and. .not. lVoluSkinSup) then
                call utmess('F', 'CHARGES_10')
            end if
            if (lVoluSkinInf .and. lVoluSkinSup) then
                ASSERT(ASTER_FALSE)
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_Mesh_module
