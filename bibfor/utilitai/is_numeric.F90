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
function is_numeric(string)
    implicit none
#include "asterf_types.h"
    character(len=*), intent(in) :: string
    aster_logical :: is_numeric
    real(kind=8) :: x
    integer(kind=8) :: e, n
    character(len=12) :: fmt
!
    n = len_trim(string)
    write (fmt, '("(f",i0,".0)")') n
    read (string, fmt, iostat=e) x
    is_numeric = e == 0
end function is_numeric
