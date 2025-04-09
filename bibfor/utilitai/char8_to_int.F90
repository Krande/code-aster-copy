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
function char8_to_int(to_convert, nommai, typent)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/is_numeric.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
!
    character(len=8), intent(in) :: to_convert
    character(len=8), optional, intent(in) :: nommai
    character(len=*), optional, intent(in) :: typent
    integer :: char8_to_int
!
    integer :: ier
    character(len=16) nomobj
    aster_logical :: lcolle
!
    lcolle = .false.
    if (present(nommai) .and. present(typent)) then
        ier = -1
        if (typent .eq. "MAILLE") then
            nomobj = nommai//".NOMMAI "
            call jeexin(nomobj, ier)
        else if (typent .eq. "NOEUD") then
            nomobj = nommai//".NOMNOE "
            call jeexin(nomobj, ier)
        else
            ASSERT(.false.)
        end if
        if (ier .eq. 0) then
            lcolle = .false.
        else
            lcolle = .true.
        end if
    end if
!
    if (.not. lcolle) then
        if (to_convert .ne. ' ') then
            if (to_convert(1:1) .eq. 'M' .or. to_convert(1:1) .eq. 'N') then
                if (.not. is_numeric(to_convert(2:8))) then
                    ASSERT(.false.)
                end if
                read (to_convert(2:8), *) char8_to_int
            else
                if (.not. is_numeric(to_convert)) then
                    ASSERT(.false.)
                end if
                read (to_convert(1:8), *) char8_to_int
            end if
        else
            char8_to_int = 0
        end if
    else
        call jenonu(jexnom(nomobj, to_convert), char8_to_int)
    end if
!
end function
