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
subroutine resuPrintIdeasElem(fileUnit, dsName, &
                              title, storeIndx, &
                              fieldType_, fieldName_, &
                              cmpUserNb_, cmpUserName_, &
                              cellUserNb_, cellUserNume_)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celcel.h"
#include "asterfort/celver.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvis.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/ircecs.h"
#include "asterfort/ircers.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/resuIdeasPermut.h"
#include "asterfort/resuSelectCmp.h"
#include "asterfort/utcmp3.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: fileUnit
    character(len=*), intent(in) :: title, dsName
    integer(kind=8), intent(in) :: storeIndx
    character(len=*), optional, intent(in) :: fieldName_, fieldType_
    integer(kind=8), optional, intent(in) :: cmpUserNb_
    character(len=8), optional, pointer :: cmpUserName_(:)
    integer(kind=8), optional, intent(in) :: cellUserNb_
    integer(kind=8), optional, pointer :: cellUserNume_(:)
!
! --------------------------------------------------------------------------------------------------
!
! Print results - IDEAS
!
! Elementary field
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
! In  cellUserNb       : number of cells require by user
! Ptr cellUserNume     : pointer to the list of index of cells require by user
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, titleLength
    character(len=8) :: meshName, quantityName
    character(len=4) :: fieldSupport
    character(len=19) :: fieldName
    character(len=16) :: fieldType
    character(len=24) :: permutJvName
    integer(kind=8), pointer :: celd(:) => null()
    character(len=24), pointer :: celk(:) => null()
    character(len=19) :: liliName
    integer(kind=8), pointer :: liel(:) => null()
    integer(kind=8), pointer :: lielLen(:) => null()
    character(len=80), pointer :: meshTitle(:) => null()
    aster_logical :: lMeshIdeas
    character(len=1) :: type
    integer(kind=8) :: fieldScalar, quantityIndx, jvVale, grelNb, iCell
    integer(kind=8) :: meshNodeNb, meshCellNb
    integer(kind=8) :: cmpCataNb, cmpListNb, cmpVariNb
    integer(kind=8), pointer :: cmpListIndx(:) => null()
    character(len=8), pointer :: cmpCataName(:) => null()
    integer(kind=8), pointer :: cmpVariIndx(:) => null()
    integer(kind=8), pointer :: permuta(:) => null()
    integer(kind=8) :: maxnod
    integer(kind=8) :: versio
    character(len=8), pointer :: cellName(:) => null()
    integer(kind=8), pointer :: cellNbNode(:) => null()
    integer(kind=8), pointer :: cellType(:) => null()
    integer(kind=8) :: cmpUserNb, cellUserNb
    character(len=8), pointer :: cmpUserName(:) => null()
    integer(kind=8), pointer :: cellUserNume(:) => null()
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
    cellUserNb = 0
    if (present(cellUserNb_)) then
        cellUserNb = cellUserNb_
        cellUserNume => cellUserNume_
    end if
!
! - For sub-points => suppress them
!
    call celver(fieldName, 'NBSPT_1', 'COOL', iret)
    if (iret .eq. 1) then
        call utmess('I', 'RESULT3_97', sk=fieldType)
        call celcel('PAS_DE_SP', fieldType, 'V', '&&IRCHML.CHAMEL2')
        fieldName = '&&IRCHML.CHAMEL2'
    end if
!
! - Access to field
!
    fieldSupport = ' '
    if (fieldName .ne. ' ') then
        call dismoi('TYPE_CHAMP', fieldName, 'CHAMP', repk=fieldSupport)
    end if
    call jeveuo(fieldName//'.CELK', 'L', vk24=celk)
    call jeveuo(fieldName//'.CELD', 'L', vi=celd)
    call jeveuo(fieldName//'.CELV', 'L', jvVale)
!
! - Physical quantity on field
!
    call jelira(fieldName//'.CELV', 'TYPE', cval=type)
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
    quantityIndx = celd(1)
    call jenuno(jexnum('&CATA.GD.NOMGD', quantityIndx), quantityName)
!
! - Access to mesh
!
    liliName = celk(1) (1:19)
    call dismoi('NOM_MAILLA', liliName, 'LIGREL', repk=meshName)
    call dismoi('NB_NO_MAILLA', meshName, 'MAILLAGE', repi=meshNodeNb)
    call dismoi('NB_MA_MAILLA', meshName, 'MAILLAGE', repi=meshCellNb)
    call jeveuo(meshName//'.TYPMAIL', 'L', vi=cellType)
!
! - Select list of components
!
    call resuSelectCmp(quantityIndx, &
                       cmpUserNb, cmpUserName, &
                       cmpCataNb, cmpCataName, &
                       cmpListNb, cmpListIndx)
    if (cmpListNb .eq. 0) then
        goto 997
    end if
!
! - Select list of components (fo VARI_R)
!
    cmpVariNb = 0
    if (cmpListNb .ne. 0) then
        if ((quantityName .eq. 'VARI_R') .and. (fieldSupport(1:2) .eq. 'EL')) then
            cmpVariNb = cmpListNb
            call utcmp3(cmpListNb, cmpUserName, cmpVariIndx)
        end if
    end if
    if (cmpListNb .eq. 0 .and. cmpUserNb .ne. 0 .and. cmpVariNb .eq. 0) then
        goto 997
    end if
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
! - SuperTab (permutation)
!
    permutJvName = '&IRCHML.PERMUTA'
    call getvis(' ', 'VERSION', scal=versio, nbret=iret)
    call resuIdeasPermut(permutJvName, versio, maxnod)
    call jeveuo(permutJvName, 'L', vi=permuta)
!
! - Get parameters of all elements in mesh
!
    AS_ALLOCATE(vk8=cellName, size=meshCellNb)
    AS_ALLOCATE(vi=cellNbNode, size=meshCellNb)
    do iCell = 1, meshCellNb
        cellName(iCell) = int_to_char8(iCell)
        call jelira(jexnum(meshName//'.CONNEX', iCell), 'LONMAX', cellNbNode(iCell))
    end do
!
! - Access to LIGREL
!
    call jeveuo(liliName//'.LIEL', 'L', vi=liel)
    call jelira(liliName//'.LIEL', 'NUTIOC', grelNb)
    call jeveuo(jexatr(liliName//'.LIEL', 'LONCUM'), 'L', vi=lielLen)
    ASSERT(grelNb .eq. celd(2))
!
! - Print elementary field
!
    if (fieldScalar .eq. 1) then
        call ircers(fileUnit, liel, grelNb, lielLen, cmpCataNb, &
                    zr(jvVale), quantityName, cmpCataName, title, cellName, &
                    fieldSupport, celd, cellNbNode, permuta, maxnod, &
                    cellType, dsName, fieldType, storeIndx, cellUserNb, &
                    cellUserNume, lMeshIdeas, cmpVariNb, cmpVariIndx, cmpListNb, &
                    cmpListIndx, cmpUserName)
    else if (fieldScalar .eq. 2) then
        call ircecs(fileUnit, liel, grelNb, lielLen, cmpCataNb, &
                    zc(jvVale), cmpCataName, title, cellName, fieldSupport, &
                    celd, cellNbNode, permuta, maxnod, cellType, &
                    dsName, fieldType, storeIndx, cellUserNb, cellUserNume, &
                    lMeshIdeas, cmpVariNb, cmpVariIndx)
    else
        ASSERT(ASTER_FALSE)
    end if
!
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
    AS_DEALLOCATE(vk8=cellName)
    AS_DEALLOCATE(vi=cellNbNode)
    AS_DEALLOCATE(vi=cmpListIndx)
    call jedetr('&&IRADHS.CODEGRA')
    call jedetr('&&IRADHS.CODEPHY')
    call jedetr('&&IRADHS.CODEPHD')
!
    call jedema()
end subroutine
