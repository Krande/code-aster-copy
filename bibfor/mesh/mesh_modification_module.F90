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
module mesh_modification_module
! ==================================================================================================
    use mesh_operators_type
    use mesh_module, only: getListOfCellGroup, getFirstNodeFromNodeGroup
    implicit none
! ==================================================================================================
    public  :: meshOperModiGetPara, meshOperModiDelPara
    public  :: meshOrieShell
    private :: meshOrieShellGetPara, meshOrieShellDelPara
! ==================================================================================================
    private
#include "asterf_types.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterc/getfac.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/ornorm.h"
#include "asterfort/orvlma.h"
#include "asterfort/utmess.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! meshOperModiGetPara
!
! Get properties from MODI_MAILLAGE command
!
! In  mesh             : name of mesh
! Out
!
! --------------------------------------------------------------------------------------------------
    subroutine meshOperModiGetPara(meshz, meshOperModiPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshz
        type(MESH_OPER_MODI_PARA), intent(out) :: meshOperModiPara
! ----- Local
        character(len=8) :: mesh
        character(len=16) :: factorKeyword
        integer(kind=8) :: nbFactorKeyword, meshDime
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        factorKeyword = "ORIE_PEAU"
        call getfac(factorKeyword, nbFactorKeyword)
        meshOperModiPara%orieSkin = nbFactorKeyword
        factorKeyword = "ORIE_NORM_COQUE"
        call getfac(factorKeyword, nbFactorKeyword)
        meshOperModiPara%orieShell = nbFactorKeyword
        factorKeyword = "ORIE_LIGNE"
        call getfac(factorKeyword, nbFactorKeyword)
        meshOperModiPara%orieLine = nbFactorKeyword
        if (meshOperModiPara%orieShell .gt. 0) then
            call meshOrieShellGetPara(meshz, meshOperModiPara)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! meshOrieShellGetPara
!
! Get parameters from MODI_MAILLAGE command
!
! In  mesh             : name of mesh
! Out meshOperModiPara : parameters from MODI_MAILLAGE command
!
! --------------------------------------------------------------------------------------------------
    subroutine meshOrieShellGetPara(meshz, meshOperModiPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshz
        type(MESH_OPER_MODI_PARA), intent(inout) :: meshOperModiPara
! ----- Local
        character(len=8) :: mesh, nodeName
        integer(kind=8) :: iFactorKeyword, nbret, nbVectCmp, nbGroupNode, nbGroupCell
        integer(kind=8) :: nodeNume, cellNume
        character(len=16), parameter :: factorKeyword = "ORIE_NORM_COQUE"
        aster_logical :: orieByVect
        real(kind=8) :: orieVect(3)
        character(len=24) :: nodeGroup, groupCellName
        integer(kind=8), pointer :: meshConnex(:) => null(), listCellNume(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ

! ----- Allocate
        allocate (meshOperModiPara%meshOperOrieShell(meshOperModiPara%orieShell))

        do iFactorKeyword = 1, meshOperModiPara%orieShell
! --------- Get VECT_NORM
            call getvr8(factorKeyword, 'VECT_NORM', iocc=iFactorKeyword, nbval=0, nbret=nbVectCmp)
            nbVectCmp = -nbVectCmp
            orieVect = 0.d0
            orieByVect = ASTER_FALSE
            if (nbVectCmp .eq. 0) then
                orieByVect = ASTER_FALSE
            else
                orieByVect = ASTER_TRUE
                call getvr8(factorKeyword, 'VECT_NORM', iocc=iFactorKeyword, nbval=nbVectCmp, &
                            vect=orieVect, nbret=nbRet)
            end if

! --------- Get GROUP_MA
            call getListOfCellGroup( &
                mesh, factorKeyword, iFactorKeyword, &
                nbGroupCell, &
                meshOperModiPara%meshOperOrieShell(iFactorKeyword)%listOfGroupOfCell)

! --------- Get GROUP_NO
            nodeGroup = " "
            nodeNume = 0
            call getvtx(factorKeyword, 'GROUP_NO', iocc=iFactorKeyword, &
                        scal=nodeGroup, nbret=nbGroupNode)
            if (nbGroupNode .ne. 0) then
                call getFirstNodeFromNodeGroup(mesh, nodeGroup, nodeName, nodeNume)
            end if
            if (nodeNume .eq. 0) then
                groupCellName = &
                    meshOperModiPara%meshOperOrieShell(iFactorKeyword)%listOfGroupOfCell(1)
                call jeveuo(jexnom(mesh//'.GROUPEMA', groupCellName), 'L', vi=listCellNume)
                cellNume = listCellNume(1)
                call jeveuo(jexnum(mesh//'.CONNEX', cellNume), 'L', vi=meshConnex)
                nodeNume = meshConnex(1)
            end if

! --------- Save parameters
            meshOperModiPara%meshOperOrieShell(iFactorKeyword)%orieByVect = orieByVect
            meshOperModiPara%meshOperOrieShell(iFactorKeyword)%orieVect = orieVect
            meshOperModiPara%meshOperOrieShell(iFactorKeyword)%nbGroupCell = nbGroupCell
            meshOperModiPara%meshOperOrieShell(iFactorKeyword)%nodeNume = nodeNume
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! meshOperModiDelPara
!
! Delete datastructure for parameters from MODI_MAILLAGE command
!
! IO  meshOperModiPara : parameters from MODI_MAILLAGE command
!
! --------------------------------------------------------------------------------------------------
    subroutine meshOperModiDelPara(meshOperModiPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_OPER_MODI_PARA), intent(inout) :: meshOperModiPara
! ----- Local
        integer(kind=8) :: iOrieShell
!   ------------------------------------------------------------------------------------------------
!
        do iOrieShell = 1, meshOperModiPara%orieShell
            call meshOrieShellDelPara(meshOperModiPara%meshOperOrieShell(iOrieShell))
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! meshOrieShellDelPara
!
! Delete datastructure for parameters from ORIE_NORM_COQUE
!
! IO  meshOperOrieShell: parameters from ORIE_NORM_COQUE
!
! --------------------------------------------------------------------------------------------------
    subroutine meshOrieShellDelPara(meshOperOrieShell)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(MESH_OPER_ORIE_SHELL), intent(inout) :: meshOperOrieShell
! ----- Local
        integer(kind=8) :: iOrieShell
!   ------------------------------------------------------------------------------------------------
!
        deallocate (meshOperOrieShell%listOfGroupOfCell)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! meshOrieShell
!
! Apply ORIE_NORM_COQUE
!
! In  mesh             : name of mesh
! IO  meshOperOrieShell: parameters from ORIE_NORM_COQUE
!
! --------------------------------------------------------------------------------------------------
    subroutine meshOrieShell(meshz, meshOperOrieShell)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshz
        type(MESH_OPER_ORIE_SHELL), intent(in) :: meshOperOrieShell
! ----- Local
        character(len=8) :: mesh
        integer(kind=8) :: iGroupCell, nbGroupCell, nbCell, nbOrie, nbOrieTotal, nodeNume, nbConex
        character(len=24) :: groupCellName
        integer(kind=8), pointer :: listCellNume(:) => null()
        aster_logical, parameter :: reorie = ASTER_TRUE
        real(kind=8) :: orieVect(3)
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        nbOrieTotal = 0

        if (meshOperOrieShell%orieByVect) then
            do iGroupCell = 1, meshOperOrieShell%nbGroupCell
                groupCellName = meshOperOrieShell%listOfGroupOfCell(iGroupCell)
                call jelira(jexnom(mesh//'.GROUPEMA', groupCellName), 'LONUTI', nbCell)
                call jeveuo(jexnom(mesh//'.GROUPEMA', groupCellName), 'L', vi=listCellNume)
                call utmess('I', "MESH3_3", sk=groupCellName, si=nbCell)
                nodeNume = meshOperOrieShell%nodeNume
                orieVect = meshOperOrieShell%orieVect
                nbOrie = 0
                call orvlma(mesh, listCellNume, nbCell, nbOrie, orieVect, nodeNume)
                nbOrieTotal = nbOrieTotal+nbOrie
                call utmess('I', "MESH3_4", si=nbOrie)
            end do
        else
            do iGroupCell = 1, meshOperOrieShell%nbGroupCell
                groupCellName = meshOperOrieShell%listOfGroupOfCell(iGroupCell)
                call jelira(jexnom(mesh//'.GROUPEMA', groupCellName), 'LONUTI', nbCell)
                call jeveuo(jexnom(mesh//'.GROUPEMA', groupCellName), 'L', vi=listCellNume)
                call utmess('I', "MESH3_3", sk=groupCellName, si=nbCell)
                nbOrie = 0
                call ornorm(mesh, listCellNume, nbCell, reorie, nbOrie, nbConex)
                if (nbConex .gt. 1) then
                    call utmess('F', 'MESH3_6')
                end if
                nbOrieTotal = nbOrieTotal+nbOrie
                call utmess('I', "MESH3_4", si=nbOrie)
            end do
        end if

        if (nbOrieTotal .ne. 0) then
            call utmess('I', "MESH3_5", si=nbOrieTotal)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module mesh_modification_module
