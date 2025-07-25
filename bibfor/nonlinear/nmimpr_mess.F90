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
subroutine nmimpr_mess(indx_mesg, vali_, valr_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/nonlinDSColumnWriteValue.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: indx_mesg
    integer(kind=8), optional, intent(in) :: vali_
    real(kind=8), optional, intent(in) :: valr_
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Print management
!
! Print UTMESS
!
! --------------------------------------------------------------------------------------------------
!
! In  indx_mesg        : index of message in MEASURE1 catalog
! In  vali             : value if integer
! In  valr             : value if real
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: time_string
!
! --------------------------------------------------------------------------------------------------
!
    if (present(valr_)) then
        call nonlinDSColumnWriteValue(0, time_=valr_, output_string_=time_string)
    end if
    if (indx_mesg .eq. 1) then
        call utmess('I', 'MEASURE1_1', sk=time_string)
    elseif (indx_mesg .eq. 2) then
        call utmess('I', 'MEASURE1_2', sk=time_string)
    elseif (indx_mesg .eq. 3) then
        call utmess('I', 'MEASURE1_3', sk=time_string)
    elseif (indx_mesg .eq. 6) then
        call utmess('I', 'MEASURE1_6', sk=time_string, si=vali_)
    elseif (indx_mesg .eq. 7) then
        call utmess('I', 'MEASURE1_7', sk=time_string, si=vali_)
    elseif (indx_mesg .eq. 8) then
        call utmess('I', 'MEASURE1_8', sk=time_string, si=vali_)
    elseif (indx_mesg .eq. 9) then
        call utmess('I', 'MEASURE1_9', sk=time_string, si=vali_)
    elseif (indx_mesg .eq. 10) then
        call utmess('I', 'MEASURE1_10', sk=time_string, si=vali_)
    elseif (indx_mesg .eq. 11) then
        call utmess('I', 'MEASURE1_11', sk=time_string)
    elseif (indx_mesg .eq. 12) then
        call utmess('I', 'MEASURE1_12', sk=time_string)
    elseif (indx_mesg .eq. 13) then
        call utmess('I', 'MEASURE1_13', sk=time_string, si=vali_)
    elseif (indx_mesg .eq. 14) then
        call utmess('I', 'MEASURE1_14', sk=time_string)
    elseif (indx_mesg .eq. 15) then
        call utmess('I', 'MEASURE1_15', sk=time_string)
    elseif (indx_mesg .eq. 16) then
        call utmess('I', 'MEASURE1_16', sk=time_string)
    elseif (indx_mesg .eq. 17) then
        call utmess('I', 'MEASURE1_17', sk=time_string)
    elseif (indx_mesg .eq. 18) then
        call utmess('I', 'MEASURE1_18', si=vali_)
    elseif (indx_mesg .eq. 19) then
        call utmess('I', 'MEASURE1_19', si=vali_)
    elseif (indx_mesg .eq. 20) then
        call utmess('I', 'MEASURE1_20', si=vali_)
    elseif (indx_mesg .eq. 21) then
        call utmess('I', 'MEASURE1_21', si=vali_)
    elseif (indx_mesg .eq. 22) then
        call utmess('I', 'MEASURE1_22', si=vali_)
    elseif (indx_mesg .eq. 23) then
        call utmess('I', 'MEASURE1_23', si=vali_)
    elseif (indx_mesg .eq. 24) then
        call utmess('I', 'MEASURE1_24', si=vali_)
    elseif (indx_mesg .eq. 25) then
        call utmess('I', 'MEASURE1_25', si=vali_)
    elseif (indx_mesg .eq. 26) then
        call utmess('I', 'MEASURE1_26', si=vali_)
    elseif (indx_mesg .eq. 27) then
        call utmess('I', 'MEASURE1_27', sk=time_string)
    elseif (indx_mesg .eq. 28) then
        call utmess('I', 'MEASURE1_28', sk=time_string)
    elseif (indx_mesg .eq. 29) then
        call utmess('I', 'MEASURE1_29', sk=time_string)
    else
        WRITE (6, *) 'IndxMesg: ', indx_mesg
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
