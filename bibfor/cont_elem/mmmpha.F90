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
subroutine mmmpha(leltf, lcont, ladhe, l_fric_no, phase)
!
    implicit none
!
#include "asterf_types.h"
!
    aster_logical, intent(in) :: leltf
    aster_logical, intent(in) :: lcont, ladhe, l_fric_no
    character(len=4), intent(out) :: phase
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Elementary computations
!
! Select phase to compute
!
! --------------------------------------------------------------------------------------------------
!
! In  leltf            : flag for friction
! In  lcont            : .true. if contact
! In  ladhe            : .true. if stick
! In  l_fric_no        : .true. if desactivation of contact during contact
! Out phase            : phase to compute
!                        'SANS' - No contact
!                        'CONT' - Contact
!                        'ADHE' - Stick
!                        'GLIS' - Slip
!                        'NCON' - Friction but not contact (!)
!
! --------------------------------------------------------------------------------------------------
!
    phase = 'SANS'
!
! - Main phase
!
    if (leltf) then
        if (lcont) then
            if (ladhe) then
                phase = 'ADHE'
            elseif (l_fric_no) then
                phase = 'NCON'
            else
                phase = 'GLIS'
            end if
        else
            phase = 'SANS'
        end if
    else
        if (lcont) then
            phase = 'CONT'
        else
            phase = 'SANS'
        end if
    end if
!
end subroutine
