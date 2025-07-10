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
function int_to_char8(to_convert, lcolle, nommai, typent)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterc/int_to_string_conversion.h"
!
    integer(kind=8), intent(in) :: to_convert
    aster_logical, optional, intent(in) :: lcolle
    character(len=8), optional, intent(in) :: nommai
    character(len=*), optional, intent(in) :: typent
    character(len=8) :: int_to_char8
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
    if (lcolle2) then
        call jenuno(jexnum(nomobj, to_convert), int_to_char8)
    else
        int_to_char8 = ' '
        call int_to_string_conversion(to_convert, int_to_char8)
    end if
!
end function
