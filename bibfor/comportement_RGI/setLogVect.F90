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

subroutine setLogVect(vect, x1, x2, x3, x4, x5, &
                      x6, x7, x8, x9, x10, x11, x12, x13, x14, &
                      x15, x16, x17)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!   copie de valeurs entières dans un tableau
!-----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
    aster_logical, intent(inout) :: vect(*)
    aster_logical, intent(in) :: x1
    aster_logical, optional, intent(in) :: x2, x3, x4, x5, x6, x7, x8
    aster_logical, optional, intent(in) :: x9, x10, x11, x12, x13, x14
    aster_logical, optional, intent(in) :: x15, x16, x17
!-----------------------------------------------------------------------
!
    vect(1) = x1
    if (present(x2)) vect(2) = x2
    if (present(x3)) vect(3) = x3
    if (present(x4)) vect(4) = x4
    if (present(x5)) vect(5) = x5
    if (present(x6)) vect(6) = x6
    if (present(x7)) vect(7) = x7
    if (present(x8)) vect(8) = x8
    if (present(x9)) vect(9) = x9
    if (present(x10)) vect(10) = x10
    if (present(x11)) vect(11) = x11
    if (present(x12)) vect(12) = x12
    if (present(x13)) vect(13) = x13
    if (present(x14)) vect(14) = x14
    if (present(x15)) vect(15) = x15
    if (present(x16)) vect(16) = x16
    if (present(x17)) vect(17) = x17
end subroutine
