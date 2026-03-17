! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine caliun(sdcontZ, meshZ, modelZ)
!
    implicit none
!
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/caraun.h"
#include "asterfort/creaun.h"
#include "asterfort/elimun.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/listun.h"
#include "asterfort/surfun.h"
#include "asterfort/wkvect.h"
#include "Contact_type.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: sdcontZ, meshZ, modelZ
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Get informations for LIAISON_UNILATER in command
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  mesh             : name of mesh
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: zoneKeyword = "ZONE"
    character(len=8) :: sdcont, mesh, model
    integer(kind=8) :: nbUnilZone, nbNodeUnil, ntCmp
    character(len=24), parameter :: nbgdcuJv = '&&CALIUN.NBGDCU', coefcuJv = '&&CARAUN.COEFCU'
    character(len=24), parameter :: compcuJv = '&&CARAUN.COMPCU', multcuJv = '&&CARAUN.MULTCU'
    character(len=24), parameter :: penacuJv = '&&CARAUN.PENACU'
    character(len=24), parameter :: noponoJv = '&&CALIUN.PONOEU', nolinoJv = '&&CALIUN.LINOEU'
    character(len=24), parameter :: lisnoeJv = '&&CALIUN.LISNOE', poinoeJv = '&&CALIUN.POINOE'
    character(len=24) :: sdunilDefi
    character(len=24) :: sdunilNdimcuJv
    integer(kind=8), pointer :: sdunilNdimcu(:) => null()
    character(len=24) :: sdcontParaciJv
    integer(kind=8), pointer :: sdcontParaci(:) => null()
    character(len=24) :: sdcontParacrJv
    real(kind=8), pointer :: sdcontParacr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    model = modelZ
    sdcont = sdcontZ
    mesh = meshZ
    sdunilDefi = sdcont(1:8)//'.UNILATE'

! - Datastructure for contact definition
    sdcontParaciJv = sdcont(1:8)//'.PARACI'
    call jeveuo(sdcontParaciJv, 'E', vi=sdcontParaci)
    sdcontParacrJv = sdcont(1:8)//'.PARACR'
    call jeveuo(sdcontParacrJv, 'E', vr=sdcontParacr)

! - Number of zones
    call getfac(zoneKeyword, nbUnilZone)
    ASSERT(nbUnilZone .ne. 0)

! - Create datastructure for general parameters
    sdunilNdimcuJv = sdunilDefi(1:16)//'.NDIMCU'
    call wkvect(sdunilNdimcuJv, 'G V I', 2, vi=sdunilNdimcu)

! - Set datastructure for general parameters
    call caraun(sdcont, zoneKeyword, nbUnilZone, &
                nbgdcuJv, coefcuJv, &
                compcuJv, multcuJv, penacuJv, ntCmp)

! - Get list of nodes
    call listun(mesh, zoneKeyword, nbUnilZone, &
                noponoJv, nolinoJv, nbNodeUnil)

! - Clean list of nodes
    call elimun(mesh, model, zoneKeyword, nbUnilZone, &
                nbgdcuJv, compcuJv, noponoJv, nolinoJv, &
                lisnoeJv, poinoeJv, &
                nbNodeUnil)

! - Get list of components
    call creaun(sdcont, mesh, model, nbUnilZone, nbNodeUnil, &
                lisnoeJv, poinoeJv, nbgdcuJv, coefcuJv, compcuJv, &
                multcuJv, penacuJv)

! - Debug
    call surfun(sdcont, mesh)

! - Clean temporary objects
    call jedetr(nolinoJv)
    call jedetr(noponoJv)
    call jedetr(lisnoeJv)
    call jedetr(poinoeJv)
    call jedetr(nbgdcuJv)
    call jedetr(coefcuJv)
    call jedetr(compcuJv)
    call jedetr(multcuJv)
    call jedetr(penacuJv)
!
    call jedema()
!
end subroutine
