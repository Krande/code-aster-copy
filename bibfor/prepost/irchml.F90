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
subroutine irchml(fileUnit, &
                  fieldTypeZ, fieldNameZ, fieldSupport, &
                  cmpUserNb, cmpUserName, &
                  cellUserNb, cellUserNume, &
                  nodeUserNb, nodeUserNume, &
                  lMeshCoor_, lmax_, lmin_, &
                  lsup_, borsup_, &
                  linf_, borinf_, &
                  realFormat_, cplxFormat_)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/celcel.h"
#include "asterfort/celces.h"
#include "asterfort/celver.h"
#include "asterfort/cesimp.h"
#include "asterfort/cncinv.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/i2trgi.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/ircecl.h"
#include "asterfort/ircerl.h"
#include "asterfort/irsspt.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/resuSelectCmp.h"
#include "asterfort/utcmp3.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: fileUnit
    character(len=*), intent(in) :: fieldNameZ, fieldTypeZ
    character(len=4), intent(in) :: fieldSupport
    integer(kind=8), intent(in) :: cmpUserNb
    character(len=8), pointer :: cmpUserName(:)
    integer(kind=8), intent(in) :: nodeUserNb
    integer(kind=8), pointer :: nodeUserNume(:)
    integer(kind=8), intent(in) :: cellUserNb
    integer(kind=8), pointer :: cellUserNume(:)
    aster_logical, optional, intent(in) :: lMeshCoor_
    aster_logical, optional, intent(in) :: lsup_, linf_, lmax_, lmin_
    real(kind=8), optional, intent(in) :: borsup_, borinf_
    character(len=*), optional, intent(in) :: realFormat_, cplxFormat_
!
! --------------------------------------------------------------------------------------------------
!
! Print results - RESULTAT
!
! Field on cells
!
! --------------------------------------------------------------------------------------------------
!
! In  fileUnit         : index of file (logical unit)
! In  fieldType        : type of field (DEPL, SIEF, EPSI, ...)
! In  fieldName        : name of field datastructure
! In  fieldSupport     : cell support of field (NOEU, ELNO, ELEM, ...)
! In  cmpUserNb        : number of components to select
! Ptr cmpUserName      : list of name of components to select
! In  nodeUserNb       : number of nodes require by user
! Ptr nodeUserNume     : list of index of nodes require by user
! In  cellUserNb       : number of cells require by user
! Ptr cellUserNume     : list of index of cells require by user
! In  lMeshCoor        : flag to print coordinates of node
! In  lmax             : flag to print maximum value on nodes
! In  lmin             : flag to print minimum value on nodes
! In  lsup             : flag if supremum exists
! In  borsup           : value of supremum
! In  linf             : flag if infinum exists
! In  borinf           : value of infinum
! In  realFormat       : format of real numbers
! In  cplxFormat       : format of complex numbers (IMAG, REAL, PHASE, MODULE or ' ')
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=8) :: meshName, quantityName
    character(len=19) :: fieldName
    character(len=19), parameter :: fieldNameS = '&&IRCHML_CES'
    character(len=16) :: fieldType
    integer(kind=8), pointer :: celd(:) => null()
    character(len=24), pointer :: celk(:) => null()
    character(len=19) :: liliName
    integer(kind=8), pointer :: liel(:) => null()
    integer(kind=8), pointer :: lielLen(:) => null()
    character(len=1) :: type
    integer(kind=8) :: fieldScalar, quantityIndx, grelNb, iCell
    character(len=24), parameter :: ncncin = '&&IRCHML.CONNECINVERSE'
    integer(kind=8) :: jdrvlc, jcncin, nbtma, iadr
    integer(kind=8) :: iNode, nodeNume, nodeNbElem, cellNumeFirst
    integer(kind=8) :: meshNodeNb, meshCellNb, meshDime
    integer(kind=8) :: cmpCataNb, cmpListNb, cmpVariNb
    integer(kind=8), pointer :: cmpListIndx(:) => null()
    character(len=8), pointer :: cmpCataName(:) => null()
    integer(kind=8), pointer :: cmpVariIndx(:) => null()
    character(len=8), pointer :: meshCellName(:) => null()
    character(len=8), pointer :: meshNodeName(:) => null()
    real(kind=8), pointer :: meshCoor(:) => null()
    integer(kind=8), pointer :: cellSelectNume(:) => null()
    integer(kind=8) :: cellSelectNb
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: connexLen(:) => null()
    aster_logical :: lMeshCoor
    aster_logical :: lsup, linf, lmax, lmin
    real(kind=8) :: borsup, borinf
    character(len=8) :: realFormat, cplxFormat
    complex(kind=8), pointer  :: valeC(:) => null()
    real(kind=8), pointer  :: valeR(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    fieldName = fieldNameZ
    fieldType = fieldTypeZ
    lMeshCoor = ASTER_FALSE
    cmpVariNb = 0

    if (present(lMeshCoor_)) then
        lMeshCoor = lMeshCoor_
    end if
    lsup = ASTER_FALSE
    if (present(lsup_)) then
        lsup = lsup_
    end if
    linf = ASTER_FALSE
    if (present(linf_)) then
        linf = linf_
    end if
    lmax = ASTER_FALSE
    if (present(lmax_)) then
        lmax = lmax_
    end if
    lmin = ASTER_FALSE
    if (present(lmin_)) then
        lmin = lmin_
    end if
    borsup = 0.d0
    if (present(borsup_)) then
        borsup = borsup_
    end if
    borinf = 0.d0
    if (present(borinf_)) then
        borinf = borinf_
    end if
    realFormat = '1PE12.5'
    if (present(realFormat_)) then
        realFormat = realFormat_
    end if
    cplxFormat = ' '
    if (present(cplxFormat_)) then
        cplxFormat = cplxFormat_
    end if
!
! - Check field "not too dynamic"
!
    call celver(fieldName, 'NBVARI_CST', 'COOL', iret)
    if (iret .eq. 1) then
        call celcel('NBVARI_CST', fieldName, 'V', '&&IRCHML.CHAMEL1')
        fieldName = '&&IRCHML.CHAMEL1'
    end if
!
! - For sub-points: special
!
    call celver(fieldName, 'NBSPT_1', 'COOL', iret)
    if (iret .eq. 1) then
        call celces(fieldName, 'V', fieldNameS)
        if (lmax .or. lmin) then
            call irsspt(fieldNameS, fileUnit, &
                        cellUserNb, cellUserNume, &
                        cmpUserNb, cmpUserName, &
                        lsup, linf, &
                        lmax, lmin, &
                        borinf, borsup)
        else
            call utmess('I', 'RESULT3_98', sk=fieldType)
            call cesimp(fieldNameS, fileUnit, cellUserNb, cellUserNume)
        end if
        call detrsd('CHAM_ELEM_S', fieldNameS)
        goto 999
    end if
!
! - Access to field
!
    call jeveuo(fieldName//'.CELK', 'L', vk24=celk)
    call jeveuo(fieldName//'.CELD', 'L', vi=celd)
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
    call dismoi('DIM_GEOM_B', meshName, 'MAILLAGE', repi=meshDime)
    call dismoi('NB_NO_MAILLA', meshName, 'MAILLAGE', repi=meshNodeNb)
    call dismoi('NB_MA_MAILLA', meshName, 'MAILLAGE', repi=meshCellNb)
    call jeveuo(meshName//'.COORDO    .VALE', 'L', vr=meshCoor)
    call jeveuo(meshName//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(meshName//'.CONNEX', 'LONCUM'), 'L', vi=connexLen)
!
! - Select list of components
!
    call resuSelectCmp(quantityIndx, &
                       cmpUserNb, cmpUserName, &
                       cmpCataNb, cmpCataName, &
                       cmpListNb, cmpListIndx)
!
! - Select list of components (for VARI_R)
!
!
    if ((quantityName .eq. 'VARI_R') .and. (fieldSupport(1:2) .eq. 'EL')) then
        if (cmpUserNb .ne. 0) then
            cmpVariNb = cmpUserNb
            AS_ALLOCATE(vi=cmpVariIndx, size=cmpVariNb)
            call utcmp3(cmpVariNb, cmpUserName, cmpVariIndx)
        end if
    end if
    if (cmpListNb .eq. 0 .and. cmpUserNb .ne. 0 .and. cmpVariNb .eq. 0) then
        goto 997
    end if
!
! - Access to LIGREL
!
    call jeveuo(liliName//'.LIEL', 'L', vi=liel)
    call jelira(liliName//'.LIEL', 'NUTIOC', grelNb)
    call jeveuo(jexatr(liliName//'.LIEL', 'LONCUM'), 'L', vi=lielLen)
    ASSERT(grelNb .eq. celd(2))
!
! - Get parameters of all elements in mesh
!
    AS_ALLOCATE(vk8=meshCellName, size=meshCellNb)
    do iCell = 1, meshCellNb
        meshCellName(iCell) = int_to_char8(iCell)
    end do
    AS_ALLOCATE(vk8=meshNodeName, size=meshNodeNb)
    do iNode = 1, meshNodeNb
        meshNodeName(iNode) = int_to_char8(iNode)
    end do
!
! - Get list of elements
!
    if (cellUserNb .eq. 0 .and. nodeUserNb .ne. 0) then
! ----- Get list of cells from list of nodes (inverse connectivity)
        call jelira(meshName//'.CONNEX', 'NMAXOC', nbtma)
        AS_ALLOCATE(vi=cellSelectNume, size=nbtma)
        call jeexin(ncncin, iret)
        if (iret .eq. 0) then
            call cncinv(meshName, [0], 0, 'V', ncncin)
        end if
        cellNumeFirst = 1
        call jeveuo(jexatr(ncncin, 'LONCUM'), 'L', jdrvlc)
        call jeveuo(jexnum(ncncin, 1), 'L', jcncin)
        do iNode = 1, nodeUserNb, 1
            nodeNume = nodeUserNume(iNode)
            nodeNbElem = zi(jdrvlc+nodeNume+1-1)-zi(jdrvlc+nodeNume-1)
            iadr = zi(jdrvlc+nodeNume-1)
            call i2trgi(cellSelectNume, zi(jcncin+iadr-1), nodeNbElem, cellNumeFirst)
        end do
        cellSelectNb = cellNumeFirst-1
    else
! ----- Get list of cells from user (trivial)
        cellSelectNb = cellUserNb
        if (cellUserNb .ne. 0) then
            AS_ALLOCATE(vi=cellSelectNume, size=cellUserNb)
            do iCell = 1, cellUserNb
                cellSelectNume(iCell) = cellUserNume(iCell)
            end do
        end if
    end if
!
! - Print elementary field
!
    if (fieldScalar .eq. 1) then
        call jeveuo(fieldName//'.CELV', 'L', vr=valeR)
        call ircerl(fileUnit, meshCellNb, liel, grelNb, lielLen, &
                    cmpCataNb, valeR, cmpCataName, meshCellName, fieldSupport, &
                    celd, connex, connexLen, meshNodeName, cmpListNb, &
                    cmpListIndx, nodeUserNb, nodeUserNume, cellSelectNb, cellSelectNume, &
                    lsup, borsup, linf, borinf, lmax, &
                    lmin, lMeshCoor, meshDime, meshCoor, liliName, &
                    realFormat, cmpVariNb, cmpVariIndx)
    else if (fieldScalar .eq. 2) then
        call jeveuo(fieldName//'.CELV', 'L', vc=valeC)
        call ircecl(fileUnit, &
                    fieldSupport, celd, realFormat, cplxFormat, &
                    nodeUserNb, nodeUserNume, &
                    cellSelectNb, cellSelectNume, &
                    meshCellNb, meshCellName, meshNodeName, &
                    lMeshCoor, meshDime, meshCoor, &
                    connex, connexLen, &
                    cmpCataNb, cmpCataName, &
                    cmpListNb, cmpListIndx, &
                    cmpVariNb, cmpVariIndx, &
                    grelNb, liel, &
                    lielLen, liliName, &
                    lmax, lmin, &
                    lsup, borsup, &
                    linf, borinf, &
                    valeC)
    else if ((fieldScalar .eq. 3) .or. (fieldScalar .eq. 4)) then
        call utmess('I', 'RESULT3_99', sk=fieldType)
        call celces(fieldName, 'V', fieldNameS)
        call cesimp(fieldNameS, fileUnit, cellSelectNb, cellSelectNume)
        call detrsd('CHAM_ELEM_S', fieldNameS)
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
    call detrsd('CHAM_ELEM', '&&IRCHML.CHAMEL1')
    call detrsd('CHAM_ELEM', '&&IRCHML.CHAMEL2')
    AS_DEALLOCATE(vi=cellSelectNume)
    AS_DEALLOCATE(vi=cmpListIndx)
    AS_DEALLOCATE(vk8=meshCellName)
    AS_DEALLOCATE(vk8=meshNodeName)
    if (cmpVariNb .ne. 0) then
        AS_DEALLOCATE(vi=cmpVariIndx)
    end if
!
    call jedema()
end subroutine
