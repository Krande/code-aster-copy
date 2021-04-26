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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine charac(load)
!
implicit none
!
#include "asterfort/adalig.h"
#include "asterfort/caddli.h"
#include "asterfort/cagene.h"
#include "asterfort/cagrou.h"
#include "asterfort/cbimpe.h"
#include "asterfort/cbvite.h"
#include "asterfort/cormgi.h"
#include "asterfort/initel.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lisnnl.h"
#include "asterfort/utmess.h"

!
character(len=8), intent(in) :: load
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads for AFFE_CHAR_ACOU_*
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: phenom = 'ACOUSTIQUE'
    character(len=4), parameter :: phenomS = 'ACOU'
    character(len=16), parameter :: command = 'AFFE_CHAR_ACOU'
    character(len=4), parameter :: valeType = 'COMP', coefType = 'REEL'
    character(len=16), parameter :: keywFactEnforceDOF = 'PRES_IMPO'
    integer :: dimeModel, iret
    character(len=8) :: mesh, model
    character(len=13) :: loadDescBase
    character(len=19) :: loadLigrel, modelLigrel
    character(len=8), pointer :: loadLigrelLgrf(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Mesh, Ligrel for model, dimension of model
    call cagene(load, command, modelLigrel, mesh, dimeModel)
    model = modelLigrel(1:8)
    if (dimeModel .gt. 3) then
        call utmess('A', 'CHARGES2_4')
    endif

! - Get Ligrel for load
    call lisnnl(phenom, load, loadDescBase)
    loadLigrel = loadDescBase//'.LIGRE'

! - Load VITE_FACE
    call cbvite(phenom, load, mesh, valeType)

! - Load IMPE_FACE
    call cbimpe(load, mesh)

! - Kinematic PRES_IMPO
    call caddli(keywFactEnforceDOF, load, mesh, modelLigrel, valeType)

! - Kinematic LIAISON_UNIF
    call cagrou(load, mesh, coefType, phenomS)

! - Update loads <LIGREL>
    call jeexin(loadLigrel//'.LGRF', iret)
    if (iret .ne. 0) then
        call adalig(loadLigrel)
        call cormgi('G', loadLigrel)
        call jeecra(loadLigrel//'.LGRF', 'DOCU', cval = phenomS)
        call initel(loadLigrel)
        call jeveuo(loadLigrel//'.LGRF', 'E', vk8 = loadLigrelLgrf)
        loadLigrelLgrf(2) = model
    endif
!
end subroutine
