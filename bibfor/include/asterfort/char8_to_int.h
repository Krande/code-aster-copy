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
#include "asterf_types.h"
!
interface
    function char8_to_int(to_convert, lcolle, nommai, typent)
        character(len=8), intent(in) :: to_convert
        aster_logical, optional, intent(in) :: lcolle
        character(len=8), optional, intent(in) :: nommai
        character(len=*), optional, intent(in) :: typent
        integer(kind=8) :: char8_to_int
    end function char8_to_int
end interface
