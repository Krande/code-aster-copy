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
subroutine irdepl(fileUnit, &
                  fieldTypeZ, fieldNameZ, &
                  cmpUserNb, cmpUserName, &
                  nodeUserNb, nodeUserNume, &
                  lMeshCoor_, lmax_, lmin_, &
                  lsup_, borsup_, &
                  linf_, borinf_, &
                  realFormat_, cplxFormat_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/ircnc8.h"
#include "asterfort/ircnrl.h"

#include "asterfort/detrsd.h"
#include "asterfort/resuSelectCmp.h"
#include "asterfort/resuSelectNode.h"
#include "asterfort/utmess.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsimp.h"
#include "asterfort/nbec.h"
!
    integer(kind=8), intent(in) :: fileUnit
    character(len=*), intent(in) :: fieldNameZ, fieldTypeZ
    integer(kind=8), intent(in) :: cmpUserNb
    character(len=8), pointer :: cmpUserName(:)
    integer(kind=8), intent(in) :: nodeUserNb
    integer(kind=8), pointer :: nodeUserNume(:)
    aster_logical, optional, intent(in) :: lMeshCoor_
    aster_logical, optional, intent(in) :: lsup_, linf_, lmax_, lmin_
    real(kind=8), optional, intent(in) :: borsup_, borinf_
    character(len=*), optional, intent(in) :: realFormat_, cplxFormat_
!
! --------------------------------------------------------------------------------------------------
!
! Print results - RESULTAT
!
! Field on nodes
!
! --------------------------------------------------------------------------------------------------
!
! In  fileUnit         : index of file (logical unit)
! In  fieldType        : type of field (DEPL, SIEF, EPSI, ...)
! In  fieldName        : name of field datastructure
! In  cmpUserNb        : number of components to select
! Ptr cmpUserName      : list of name of components to select
! In  nodeUserNb       : number of nodes require by user
! Ptr nodeUserNume     : list of index of nodes require by user
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
    integer(kind=8) :: cmpCataNb, cmpListNb
    character(len=8) :: meshName
    character(len=19) :: fieldName
    character(len=16) :: fieldType
    character(len=19), parameter :: fieldNameS = '&&IRDEPL_CES'
    character(len=24) :: profName
    character(len=1) :: type
    integer(kind=8) :: fieldScalar, quantityIndx, nec, liliMesh
    integer(kind=8) :: meshDime, meshNodeNb
    integer(kind=8), pointer :: cmpListIndx(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    integer(kind=8), pointer :: prno(:) => null()
    real(kind=8), pointer :: meshCoor(:) => null()
    character(len=8), pointer :: cmpCataName(:) => null()
    integer(kind=8) :: nodeNb
    character(len=8), pointer :: nodeListName(:) => null()
    integer(kind=8), pointer :: nodeListNume(:) => null()
    integer(kind=8), pointer :: codeInte(:) => null()
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
! - Get properties of field
!
    call jelira(fieldName//'.VALE', 'TYPE', cval=type)
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
    call dismoi('NUM_GD', fieldName, 'CHAM_NO', repi=quantityIndx)
!
! - "coded" integers
!
    nec = nbec(quantityIndx)
    AS_ALLOCATE(vi=codeInte, size=nec)
!
! - Access to mesh
!
    call dismoi('NOM_MAILLA', fieldName, 'CHAM_NO', repk=meshName)
    call dismoi('DIM_GEOM_B', meshName, 'MAILLAGE', repi=meshDime)
    call jeveuo(meshName//'.COORDO    .VALE', 'L', vr=meshCoor)
    call dismoi('NB_NO_MAILLA', meshName, 'MAILLAGE', repi=meshNodeNb)
!
! - Access to profile of numbering
!
    call dismoi('NUME_EQUA', fieldName, 'CHAM_NO', repk=profName)
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
    AS_ALLOCATE(vk8=nodeListName, size=meshNodeNb)
    AS_ALLOCATE(vi=nodeListNume, size=meshNodeNb)
    call resuSelectNode(meshName, meshNodeNb, &
                        nodeUserNb, nodeUserNume, &
                        nodeListName, nodeListNume, &
                        nodeNb)
!
! - Print nodal field
!
    if (fieldScalar .eq. 1) then
        call jeveuo(fieldName//'.VALE', 'L', vr=valeR)
        call ircnrl(fileUnit, nodeNb, prno, nueq, nec, &
                    codeInte, cmpCataNb, valeR, cmpCataName, nodeListName, &
                    lMeshCoor, meshDime, meshCoor, nodeListNume, cmpListNb, &
                    cmpListIndx, lsup, borsup, linf, borinf, &
                    lmax, lmin, realFormat)
    else if (fieldScalar .eq. 2) then
        call jeveuo(fieldName//'.VALE', 'L', vc=valeC)
        call ircnc8(fileUnit, realFormat, cplxFormat, &
                    nodeNb, nodeListNume, nodeListName, &
                    lMeshCoor, meshDime, meshCoor, &
                    cmpCataNb, cmpCataName, &
                    cmpListNb, cmpListIndx, &
                    nec, nueq, &
                    prno, codeInte, &
                    lmax, lmin, &
                    lsup, borsup, &
                    linf, borinf, &
                    valeC)
    else if (fieldScalar .eq. 3 .or. fieldScalar .eq. 4) then
        call utmess('I', 'RESULT3_99', sk=fieldType)
        call cnocns(fieldName, 'V', fieldNameS)
        call cnsimp(fieldNameS, fileUnit)
        call detrsd('CHAM_NO_S', fieldNameS)
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
    AS_DEALLOCATE(vi=cmpListIndx)
    AS_DEALLOCATE(vi=codeInte)
    AS_DEALLOCATE(vk8=nodeListName)
    AS_DEALLOCATE(vi=nodeListNume)
!
    call jedema()
end subroutine
