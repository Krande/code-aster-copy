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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine hasBehaviourFeature(model, caraElem, compor, feature, lFlag)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesred.h"
#include "asterfort/cesvar.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jeveuo.h"
!
    character(len=8), intent(in) :: model, caraElem
    character(len=19), intent(in) :: compor
    character(len=16), intent(in) :: feature
    aster_logical, intent(out) :: lFlag
!
! --------------------------------------------------------------------------------------------------
!
! Constitutive laws
!
! Detecting the presence of certain features in behaviour
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  caraElem         : name of the map for elementary characteristics
! In  compor           : name of the map for behaviour
! In  feature          : feature to detect
! Out lFlag            : flag for presence of feature
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=24) :: modelLigrel, cmpName
    character(len=16) :: valueName
    character(len=19), parameter :: comporS = '&&NMCPQU.COTO'
    character(len=19), parameter :: comporR = '&&NMCPQU.COPM'
    integer(kind=8) :: nbCell, iCell
    integer(kind=8) :: jdecal, jcesd, jcesl
    character(len=16), pointer :: cesv(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    lFlag = ASTER_FALSE

! - Create a "simplified" field to extend behaviour field (DCEL_I)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call exisd('CHAM_ELEM_S', compor, iret)
    if (iret .eq. 0) then
        call cesvar(caraElem, compor, modelLigrel, compor)
    end if

! - Name of component and value in field for feature
    cmpName = " "
    valueName = " "
    if (feature .eq. "Annealing") then
        cmpName = "POSTINCR"
        valueName = "REST_ECRO"
    else
        ASSERT(ASTER_FALSE)
    end if

! - Transform map in SimpleFieldOnCells
    call carces(compor, 'ELEM', ' ', 'V', comporS, 'A', iret)

! - Reduced field on component
    call cesred(comporS, 0, [0], 1, cmpName, 'V', comporR)

! - Access to reduced simpleFieldOnCells
    call jeveuo(comporR(1:19)//'.CESD', 'L', jcesd)
    call jeveuo(comporR(1:19)//'.CESL', 'L', jcesl)
    call jeveuo(comporR(1:19)//'.CESV', 'L', vk16=cesv)
    nbCell = zi(jcesd-1+1)

! - Look for value in components
    do iCell = 1, nbCell
        call cesexi('C', jcesd, jcesl, iCell, 1, 1, 1, jdecal)
        ASSERT(jdecal .gt. 0)
        if (cesv(jdecal) .eq. valueName) then
            lFlag = ASTER_TRUE
            goto 99
        end if
    end do
!
99  continue
!
end subroutine
