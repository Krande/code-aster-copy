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
subroutine resuPrintIdeasNode(fileUnit, dsName, &
                              title, storeIndx, &
                              fieldType_, fieldName_, &
                              cmpUserNb_, cmpUserName_, &
                              nodeUserNb_, nodeUserNume_)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/irdesc.h"
#include "asterfort/irdesr.h"

#include "asterfort/resuSelectCmp.h"
#include "asterfort/resuSelectNode.h"
#include "asterfort/utmess.h"
#include "asterfort/nbec.h"
!
    integer(kind=8), intent(in) :: fileUnit
    character(len=*), intent(in) :: title, dsName
    integer(kind=8), intent(in) :: storeIndx
    character(len=*), optional, intent(in) :: fieldName_, fieldType_
    integer(kind=8), optional, intent(in) :: cmpUserNb_
    character(len=8), optional, pointer :: cmpUserName_(:)
    integer(kind=8), optional, intent(in) :: nodeUserNb_
    integer(kind=8), optional, pointer :: nodeUserNume_(:)
!
! --------------------------------------------------------------------------------------------------
!
! Print results - IDEAS
!
! Nodal field
!
! --------------------------------------------------------------------------------------------------
!
! In  fileUnit         : index of file (logical unit)
! In  dsName           : name of datastructure (result or field)
! In  title            : title of result
! In  fieldType        : type of field (DEPL, SIEF, EPSI, ...)
! In  fieldName        : name of field datastructure
! In  storeIndx        : index of slot to save field
! In  cmpUserNb        : number of components to select
! Ptr cmpUserName      : pointer to the names of components to select
! In  nodeUserNb       : number of nodes require by user
! Ptr nodeUserNume     : pointer to the list of index of nodes require by user
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, titleLength
    integer(kind=8) :: cmpCataNb, cmpListNb
    character(len=8) :: meshName
    character(len=19) :: fieldName
    character(len=16) :: fieldType
    character(len=24) :: profName
    character(len=80), pointer :: meshTitle(:) => null()
    aster_logical :: lMeshIdeas
    character(len=1) :: type
    integer(kind=8) :: fieldScalar, quantityIndx, nec, liliMesh, jvVale
    integer(kind=8) :: meshDime, meshNodeNb
    integer(kind=8), pointer :: cmpListIndx(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    integer(kind=8), pointer :: prno(:) => null()
    character(len=8), pointer :: cmpCataName(:) => null()
    integer(kind=8) :: nodeNb
    character(len=8), pointer :: nodeName(:) => null()
    integer(kind=8), pointer :: nodeNume(:) => null()
    integer(kind=8), pointer :: codeInte(:) => null()
    integer(kind=8) :: cmpUserNb, nodeUserNb
    character(len=8), pointer :: cmpUserName(:) => null()
    integer(kind=8), pointer :: nodeUserNume(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    fieldName = ' '
    if (present(fieldName_)) then
        fieldName = fieldName_
    end if
    fieldType = ' '
    if (present(fieldType_)) then
        fieldType = fieldType_
    end if
    cmpUserNb = 0
    if (present(cmpUserNb_)) then
        cmpUserNb = cmpUserNb_
        cmpUserName => cmpUserName_
    end if
    nodeUserNb = 0
    if (present(nodeUserNb_)) then
        nodeUserNb = nodeUserNb_
        nodeUserNume => nodeUserNume_
    end if
!
! - Get properties of field
!
    call jelira(fieldName//'.VALE', 'TYPE', cval=type)
    call jeveuo(fieldName//'.VALE', 'L', jvVale)
    if (type(1:1) .eq. 'R') then
        fieldScalar = 1
    else if (type(1:1) .eq. 'C') then
        fieldScalar = 2
    else if (type(1:1) .eq. 'I') then
        fieldScalar = 3
    else if (type(1:1) .eq. 'K') then
        fieldScalar = 4
    else
        ASSERT(ASTER_FALSE)
    end if
    call dismoi("NUM_GD", fieldName, "CHAM_NO", repi=quantityIndx)
    call dismoi("NOM_MAILLA", fieldName, "CHAM_NO", repk=meshName)
    call dismoi("NUME_EQUA", fieldName, "CHAM_NO", repk=profName)
!
! - "coded" integers
!
    nec = nbec(quantityIndx)
    AS_ALLOCATE(vi=codeInte, size=nec)
!
! - Access to mesh
!
    call dismoi('DIM_GEOM_B', meshName, 'MAILLAGE', repi=meshDime)
    call dismoi('NB_NO_MAILLA', meshName, 'MAILLAGE', repi=meshNodeNb)
!
! - Access to profile of numbering
!
    call jeveuo(profName(1:19)//'.NUEQ', 'L', vi=nueq)
    call jenonu(jexnom(profName(1:19)//'.LILI', '&MAILLA'), liliMesh)
    call jeveuo(jexnum(profName(1:19)//'.PRNO', liliMesh), 'L', vi=prno)
!
! - Select list of components
!
    call resuSelectCmp(quantityIndx, &
                       cmpUserNb, cmpUserName, &
                       cmpCataNb, cmpCataName, &
                       cmpListNb, cmpListIndx)
    if (cmpListNb .eq. 0 .and. cmpUserNb .ne. 0) then
        goto 997
    end if
!
! - Select list of nodes
!
    AS_ALLOCATE(vk8=nodeName, size=meshNodeNb)
    AS_ALLOCATE(vi=nodeNume, size=meshNodeNb)
    call resuSelectNode(meshName, meshNodeNb, &
                        nodeUserNb, nodeUserNume, &
                        nodeName, nodeNume, &
                        nodeNb)
!
! - Is an IDEAS mesh ?
!
    lMeshIdeas = ASTER_FALSE
    call jeexin(meshName//'           .TITR', iret)
    if (iret .ne. 0) then
        call jeveuo(meshName//'           .TITR', 'L', vk80=meshTitle)
        call jelira(meshName//'           .TITR', 'LONMAX', titleLength)
        if (titleLength .ge. 1) then
            if (meshTitle(1) (10:31) .eq. 'AUTEUR=INTERFACE_IDEAS') then
                lMeshIdeas = ASTER_TRUE
            end if
        end if
    end if
!
! - Print nodal field
!
    if (fieldScalar .eq. 1) then
        call irdesr(fileUnit, nodeNb, prno, nueq, nec, &
                    codeInte, cmpCataNb, zr(jvVale), cmpCataName, title, &
                    nodeName, dsName, fieldType, storeIndx, nodeNume, &
                    lMeshIdeas, cmpListNb, cmpListIndx, cmpUserName)
    else if (fieldScalar .eq. 2) then
        call irdesc(fileUnit, nodeNb, prno, nueq, nec, &
                    codeInte, cmpCataNb, zc(jvVale), cmpCataName, title, &
                    nodeName, dsName, fieldType, storeIndx, nodeNume, &
                    lMeshIdeas)
    end if
    goto 999
!
997 continue
!
    call utmess('A', 'RESULT3_40', sk=fieldType)
!
999 continue
!
! - Clean
!
    AS_DEALLOCATE(vi=codeInte)
    AS_DEALLOCATE(vk8=nodeName)
    AS_DEALLOCATE(vi=nodeNume)
    AS_DEALLOCATE(vi=cmpListIndx)
!
    call jedema()
end subroutine
