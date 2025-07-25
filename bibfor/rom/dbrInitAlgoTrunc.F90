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
subroutine dbrInitAlgoTrunc(paraTrunc)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/dismoi.h"
#include "asterfort/gnomsd.h"
#include "asterfort/infniv.h"
#include "asterfort/modelNodeEF.h"
#include "asterfort/numero.h"
#include "asterfort/romFieldNodesAreDefined.h"
#include "asterfort/addModelLigrel.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_ParaDBR_Trunc), intent(inout) :: paraTrunc
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_BASE_REDUITE
!
! Initializations for algorith - For truncation
!
! --------------------------------------------------------------------------------------------------
!
! IO  paraTrunc         : datastructure for parameters (truncation)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbLigr
    character(len=24), pointer :: listLigr(:) => null()
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbEquaRom, nbNodeRom
    character(len=8) :: modelRom, modelDom
    character(len=24) :: numeRom, numeDom, noojb
    integer(kind=8), pointer :: numeNodeRom(:) => null()
    type(ROM_DS_Field) :: mode
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM18_32')
    end if
!
! - Get parameters
!
    noojb = '12345678.00000.NUME.PRNO'
    mode = paraTrunc%baseInit%mode
    modelRom = paraTrunc%modelRom
    modelDom = paraTrunc%baseInit%mode%model

! - Create name of object for numbering for ROM
    if (niv .ge. 2) then
        call utmess('I', 'ROM18_33')
    end if
    numeRom = '12345678.NUMED'
    call gnomsd(' ', noojb, 10, 14)
    numeRom = noojb(1:14)

! - Add LIGREL from model
    nbLigr = 0
    call addModelLigrel(modelRom, nbLigr, listLigr)

! - Create numbering for ROM
    call numero(numeRom, 'VV', &
                nbLigr, listLigr)
    AS_DEALLOCATE(vk24=listLigr)

! - Create name of object for numbering for DOM
    if (niv .ge. 2) then
        call utmess('I', 'ROM18_34')
    end if
    numeDom = '12345678.NUMED'
    call gnomsd(' ', noojb, 10, 14)
    numeDom = noojb(1:14)

! - Add LIGREL from model
    nbLigr = 0
    call addModelLigrel(modelDom, nbLigr, listLigr)

! - Create numbering for DOM
    call numero(numeDom, 'VV', &
                nbLigr, listLigr)
    AS_DEALLOCATE(vk24=listLigr)

! - Extract list of nodes on reduced model
    call modelNodeEF(modelRom, nbNodeRom, numeNodeRom)

! - Prepare the list of equations from list of nodes
    call romFieldNodesAreDefined(mode, paraTrunc%equaRom, numeDom, &
                                 nbNode_=nbNodeRom, &
                                 listNode_=numeNodeRom)

! - Save parameters
    call dismoi('NB_EQUA', numeRom, 'NUME_DDL', repi=nbEquaRom)
    paraTrunc%nbEquaRom = nbEquaRom
!
end subroutine
