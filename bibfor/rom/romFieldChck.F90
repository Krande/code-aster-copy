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
subroutine romFieldChck(field, fieldName_)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
    type(ROM_DS_Field), intent(in) :: field
    character(len=*), optional, intent(in) :: fieldName_
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Field management
!
! Check components
!
! --------------------------------------------------------------------------------------------------
!
! In  field            : field
! In  fieldName        : name of field where empiric modes have been constructed
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: fieldName
    integer(kind=8) :: nbCmpChck, nbCmpName
    integer(kind=8) :: iCmpName, iCmpChck, cmpIndx
    character(len=8) :: chckCmpName(6), cmpName
!
! --------------------------------------------------------------------------------------------------
!
    fieldName = ' '
    nbCmpName = field%nbCmpName
    if (present(fieldName_)) then
        fieldName = fieldName_
    else
        fieldName = field%fieldName
    end if
!
! - List of components authorized in field
!
    if (fieldName .eq. 'TEMP') then
        nbCmpChck = 1
        chckCmpName(1) = 'TEMP'
    elseif (fieldName .eq. 'DEPL') then
        nbCmpChck = 3
        chckCmpName(1) = 'DX'
        chckCmpName(2) = 'DY'
        chckCmpName(3) = 'DZ'
    elseif (fieldName .eq. 'FLUX_NOEU') then
        nbCmpChck = 3
        chckCmpName(1) = 'FLUX'
        chckCmpName(2) = 'FLUY'
        chckCmpName(3) = 'FLUZ'
    elseif (fieldName .eq. 'SIEF_NOEU') then
        nbCmpChck = 6
        chckCmpName(1) = 'SIXX'
        chckCmpName(2) = 'SIYY'
        chckCmpName(3) = 'SIZZ'
        chckCmpName(4) = 'SIXZ'
        chckCmpName(5) = 'SIYZ'
        chckCmpName(6) = 'SIXY'
    elseif (fieldName .eq. 'SIEF_ELGA') then
        nbCmpChck = 6
        chckCmpName(1) = 'SIXX'
        chckCmpName(2) = 'SIYY'
        chckCmpName(3) = 'SIZZ'
        chckCmpName(4) = 'SIXZ'
        chckCmpName(5) = 'SIYZ'
        chckCmpName(6) = 'SIXY'
    elseif (fieldName .eq. 'UPPHI_2D') then
        nbCmpChck = 4
        chckCmpName(1) = 'DX'
        chckCmpName(2) = 'DY'
        chckCmpName(3) = 'PRES'
        chckCmpName(4) = 'PHI'
    elseif (fieldName .eq. 'UPPHI_3D') then
        nbCmpChck = 5
        chckCmpName(1) = 'DX'
        chckCmpName(2) = 'DY'
        chckCmpName(3) = 'DZ'
        chckCmpName(4) = 'PRES'
        chckCmpName(5) = 'PHI'
    elseif (fieldName .eq. 'VARI_ELGA') then
        nbCmpChck = 0
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Required components
!
    do iCmpChck = 1, nbCmpChck
        cmpName = chckCmpName(iCmpChck)
        cmpIndx = indik8(field%listCmpName, chckCmpName(iCmpChck), 1, nbCmpName)
        if (cmpIndx .eq. 0) then
            call utmess('F', 'ROM11_25', sk=cmpName)
        end if
    end do
!
! - Forbidden components
!
    if (nbCmpChck .gt. 0) then
        do iCmpName = 1, nbCmpName
            cmpName = field%listCmpName(iCmpName)
            cmpIndx = indik8(chckCmpName, cmpName, 1, nbCmpChck)
            if (cmpIndx .eq. 0) then
                call utmess('F', 'ROM11_23', sk=cmpName)
            end if
        end do
    end if
!
end subroutine
