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
!
#include "asterf_types.h"
!
interface
    function meiden(scal, ncmp, i1, i3, nec,&
                    i2, i4)
        character(len=4) :: scal
        integer(kind=8) :: ncmp
        integer(kind=8) :: i1
        integer(kind=8) :: i3
        integer(kind=8) :: nec
        integer(kind=8) :: i2
        integer(kind=8) :: i4
        aster_logical :: meiden
    end function meiden
end interface
