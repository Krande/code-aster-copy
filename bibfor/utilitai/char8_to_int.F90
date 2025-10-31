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
function char8_to_int(to_convert, lcolle, nommai, typent)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/is_numeric.h"
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
!
    character(len=8), intent(in) :: to_convert
    aster_logical, optional, intent(in) :: lcolle
    character(len=8), optional, intent(in) :: nommai
    character(len=*), optional, intent(in) :: typent
    integer(kind=8) :: char8_to_int
!
    character(len=16) nomobj
    aster_logical :: lcolle2
!
    if (present(lcolle)) then
        lcolle2 = lcolle
    else
        lcolle2 = .false.
    end if
    if (present(nommai) .and. present(typent)) then
        if (typent .eq. "MAILLE") then
            nomobj = nommai//".NOMMAI "
        else if (typent .eq. "NOEUD") then
            nomobj = nommai//".NOMNOE "
        else
            ASSERT(.false.)
        end if
    end if
!
    if (.not. lcolle2) then
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
