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
subroutine carc_init(mesh, carcri, base)
!
    implicit none

#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: carcri
    character(len=1), intent(in) :: base
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Initialization of map
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : mesh
! In  carcri           : map for parameters for integration of constitutive law
! In  base             : Jeveux base
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: physQuantityName = 'CARCRI'
    integer(kind=8) :: physQuantityNume
    integer(kind=8) :: nbCmp, iCmp
    real(kind=8), pointer :: carcriValv(:) => null()
    character(len=8), pointer :: cataNomcmp(:) => null()
    character(len=8), pointer :: carcriNcmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!

! - Read catalog
    call jenonu(jexnom('&CATA.GD.NOMGD', physQuantityName), physQuantityNume)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', physQuantityNume), 'L', vk8=cataNomcmp)
    call jelira(jexnum('&CATA.GD.NOMCMP', physQuantityNume), 'LONMAX', nbCmp)
    ASSERT(nbCmp .eq. CARCRI_SIZE)

! - Allocate <CARTE>
    call detrsd("CARTE", carcri)
    call alcart(base, carcri, mesh, physQuantityName)

! - Acces to <CARTE>
    call jeveuo(carcri(1:19)//'.NCMP', 'E', vk8=carcriNcmp)
    call jeveuo(carcri(1:19)//'.VALV', 'E', vr=carcriValv)

! - Init <CARTE>
    do iCmp = 1, nbCmp
        carcriNcmp(iCmp) = cataNomcmp(iCmp)
        carcriValv(iCmp) = 0.d0
    end do

! - Default values
    carcriValv(1) = 10
    carcriValv(2) = 0
    carcriValv(3) = 1.d-6
    carcriValv(4) = 1.d0
    carcriValv(4) = 1.d0

! - Allocate on all mesh
    call nocart(carcri, 1, nbCmp)
!
end subroutine
