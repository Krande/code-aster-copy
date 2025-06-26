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
subroutine checkSuperElement(optionZ, modelZ)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: optionZ, modelZ
!
! --------------------------------------------------------------------------------------------------
!
! Check if super-elements have been computed
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : option to compute
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: option
    character(len=8) :: model, mesh, superElementName
    integer(kind=8) :: iCell, iexi
    integer(kind=8) :: nbSuperElement, nbSuperCell
    integer(kind=8), pointer :: sssa(:) => null()
    character(len=8), pointer :: nomacr(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    option = optionZ
    model = modelZ
    ASSERT(model(1:1) .ne. ' ')
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nbSuperElement)
    if (nbSuperElement .gt. 0) then
        call dismoi('NB_SM_MAILLA', model, 'MODELE', repi=nbSuperCell)
        call jeveuo(mesh//'.NOMACR', 'L', vk8=nomacr)
        call jeveuo(model//'.MODELE    .SSSA', 'L', vi=sssa)
        if (option .eq. 'MASS_MECA') then
            do iCell = 1, nbSuperCell
                superElementName = nomacr(iCell)
                if (sssa(iCell) .eq. 1) then
                    call jeexin(superElementName//'.MAEL_MASS_VALE', iexi)
                    if (iexi .eq. 0) then
                        call utmess('E', 'SUPERELEMENT1_1')
                    end if
                end if
            end do
        else if (option .eq. 'RIGI_MECA') then
            do iCell = 1, nbSuperCell
                superElementName = nomacr(iCell)
                if (sssa(iCell) .eq. 1) then
                    call jeexin(superElementName//'.MAEL_RAID_VALE', iexi)
                    if (iexi .eq. 0) then
                        call utmess('E', 'SUPERELEMENT1_2')
                    end if
                end if
            end do
        else if (option .eq. 'AMOR_MECA') then
            do iCell = 1, nbSuperCell
                superElementName = nomacr(iCell)
                if (sssa(iCell) .eq. 1) then
                    call jeexin(superElementName//'.MAEL_AMOR_VALE', iexi)
                    if (iexi .eq. 0) then
                        call utmess('F', 'SUPERELEMENT1_3')
                    end if
                end if
            end do
        end if
    end if
!
end subroutine
