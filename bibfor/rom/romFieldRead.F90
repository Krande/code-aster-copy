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
subroutine romFieldRead(operation, field, fieldObject, &
                        fieldVale_, resultName_, numeStore_)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jelibe.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: operation
    type(ROM_DS_Field), intent(in) :: field
    character(len=24), intent(inout) :: fieldObject
    character(len=8), optional, intent(in) :: resultName_
    integer(kind=8), optional, intent(in) :: numeStore_
    real(kind=8), optional, pointer :: fieldVale_(:)
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Field management
!
! Read field from results datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  operation        : operation
!                        'Read'      to read field
!                        'Free'      to free memory from field
! In  field            : structure of field to read
! In  fieldObject      : name of JEVEUX object of field
! In  resultName       : name of results datastructure to read
! In  numeStore        : index to read field in results datastructure
! Ptr fieldVale        : pointer for values of field
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret, iEqua, nbEqua
    character(len=24) :: fieldName
    character(len=4) :: fieldSupp
    aster_logical :: lFilter
!
! --------------------------------------------------------------------------------------------------
!
    fieldName = field%fieldName
    fieldSupp = field%fieldSupp
    lFilter = field%lFilter
    nbEqua = field%nbEqua
!
    if (operation .eq. 'Read') then
        call rsexch(' ', resultName_, fieldName, &
                    numeStore_, fieldObject, iret)
        if (iret .ne. 0) then
            call utmess('F', 'ROM11_11', sk=fieldName, si=numeStore_)
        end if
        if (fieldSupp == 'NOEU') then
            call jeveuo(fieldObject(1:19)//'.VALE', 'L', vr=fieldVale_)
        else if (fieldSupp == 'ELGA') then
            call jeveuo(fieldObject(1:19)//'.CELV', 'L', vr=fieldVale_)
        else
            ASSERT(ASTER_FALSE)
        end if
        if (lFilter) then
            do iEqua = 1, nbEqua
                fieldVale_(iEqua) = fieldVale_(iEqua)*field%equaFilter(iEqua)
            end do
        end if

    elseif (operation .eq. 'Free') then
        if (fieldSupp == 'NOEU') then
            call jelibe(fieldObject(1:19)//'.VALE')
        elseif (fieldSupp == 'ELGA') then
            call jelibe(fieldObject(1:19)//'.CELV')
        else
            ASSERT(ASTER_FALSE)
        end if

    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
