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
    subroutine getLogVect(vect, x1, x2, x3, x4, x5,&
                      x6, x7, x8, x9, x10, x11, x12, x13, x14,&
                      x15, x16, x17)
        aster_logical, intent(in) :: vect(*)
        aster_logical, intent(out) :: x1
        aster_logical, optional, intent(out) :: x2
        aster_logical, optional, intent(out) :: x3
        aster_logical, optional, intent(out) :: x4
        aster_logical, optional, intent(out) :: x5
        aster_logical, optional, intent(out) :: x6
        aster_logical, optional, intent(out) :: x7
        aster_logical, optional, intent(out) :: x8
        aster_logical, optional, intent(out) :: x9
        aster_logical, optional, intent(out) :: x10
        aster_logical, optional, intent(out) :: x11
        aster_logical, optional, intent(out) :: x12
        aster_logical, optional, intent(out) :: x13
        aster_logical, optional, intent(out) :: x14
        aster_logical, optional, intent(out) :: x15
        aster_logical, optional, intent(out) :: x16
        aster_logical, optional, intent(out) :: x17
    end subroutine getLogVect
end interface
