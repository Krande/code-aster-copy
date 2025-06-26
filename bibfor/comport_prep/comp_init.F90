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
subroutine comp_init(mesh, compor, base, nbCmp_)
!
    implicit none
!
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: compor
    character(len=1), intent(in) :: base
    integer(kind=8), optional, intent(out) :: nbCmp_
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviours  (all)
!
! Initialization of COMPOR <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : mesh
! In  compor           : map for parameters of constitutive laws
! In  base             : permanent or temporary allocation
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: physQuantityName = 'COMPOR'
    integer(kind=8) :: physQuantityNume
    integer(kind=8) :: nbCmp, iCmp
    character(len=16), pointer :: comporValv(:) => null()
    character(len=8), pointer :: cataNomcmp(:) => null()
    character(len=8), pointer :: comporNcmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Read catalog
    call jenonu(jexnom('&CATA.GD.NOMGD', physQuantityName), physQuantityNume)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', physQuantityNume), 'L', vk8=cataNomcmp)
    call jelira(jexnum('&CATA.GD.NOMCMP', physQuantityNume), 'LONMAX', nbCmp)
    ASSERT(nbCmp .le. COMPOR_SIZE)

! - Allocate <CARTE>
    call detrsd("CARTE", compor)
    call alcart(base, compor, mesh, physQuantityName)

! - Acces to <CARTE>
    call jeveuo(compor(1:19)//'.NCMP', 'E', vk8=comporNcmp)
    call jeveuo(compor(1:19)//'.VALV', 'E', vk16=comporValv)

! - Init <CARTE>
    do iCmp = 1, nbCmp
        comporNcmp(iCmp) = cataNomcmp(iCmp)
        comporValv(iCmp) = ' '
    end do

! - Output
    if (present(nbCmp_)) then
        nbCmp_ = nbCmp
    end if
!
end subroutine
