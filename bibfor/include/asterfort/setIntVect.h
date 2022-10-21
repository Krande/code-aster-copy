! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

! --------------------------------------------------------------------
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANvectILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
    subroutine setIntVect(vect, x1, x2, x3, x4, x5,&
                      x6, x7, x8, x9, x10, x11, x12, x13, x14,&
                      x15, x16, x17)
        integer, intent(inout) :: vect(*)
        integer, intent(in) :: x1
        integer, optional, intent(in) :: x2
        integer, optional, intent(in) :: x3
        integer, optional, intent(in) :: x4
        integer, optional, intent(in) :: x5
        integer, optional, intent(in) :: x6
        integer, optional, intent(in) :: x7
        integer, optional, intent(in) :: x8
        integer, optional, intent(in) :: x9
        integer, optional, intent(in) :: x10
        integer, optional, intent(in) :: x11
        integer, optional, intent(in) :: x12
        integer, optional, intent(in) :: x13
        integer, optional, intent(in) :: x14
        integer, optional, intent(in) :: x15
        integer, optional, intent(in) :: x16
        integer, optional, intent(in) :: x17
    end subroutine setIntVect
end interface
