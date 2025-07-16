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

subroutine setValVect(vect, x1, x2, x3, x4, x5, &
                      x6, x7, x8, x9, x10, x11, x12, x13, x14, &
                      x15, x16, x17, ind1, vectInd)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!   copie de valeur dans un tableau
!   soit on donne le premier indice ind1 et les autres suivent
!   soit on donne le vecteur d'indice vectInd
!-----------------------------------------------------------------------
    implicit none
#include "asterfort/assert.h"
    real(kind=8), intent(inout) :: vect(*)
    real(kind=8), intent(in) :: x1
    real(kind=8), optional, intent(in) :: x2, x3, x4, x5, x6, x7, x8
    real(kind=8), optional, intent(in) :: x9, x10, x11, x12, x13, x14
    real(kind=8), optional, intent(in) :: x15, x16, x17
    integer(kind=8), optional, intent(in) :: ind1, vectInd(:)
!-----------------------------------------------------------------------
!
    ASSERT(present(ind1) .or. present(vectInd))
!
    if (present(ind1)) then
        ASSERT(ind1 .ge. 1)
        vect(ind1) = x1
        if (present(x2)) vect(ind1-1+2) = x2
        if (present(x3)) vect(ind1-1+3) = x3
        if (present(x4)) vect(ind1-1+4) = x4
        if (present(x5)) vect(ind1-1+5) = x5
        if (present(x6)) vect(ind1-1+6) = x6
        if (present(x7)) vect(ind1-1+7) = x7
        if (present(x8)) vect(ind1-1+8) = x8
        if (present(x9)) vect(ind1-1+9) = x9
        if (present(x10)) vect(ind1-1+10) = x10
        if (present(x11)) vect(ind1-1+11) = x11
        if (present(x12)) vect(ind1-1+12) = x12
        if (present(x13)) vect(ind1-1+13) = x13
        if (present(x14)) vect(ind1-1+14) = x14
        if (present(x15)) vect(ind1-1+15) = x15
        if (present(x16)) vect(ind1-1+16) = x16
        if (present(x17)) vect(ind1-1+17) = x17
    elseif (present(vectInd)) then
        vect(vectInd(1)) = x1
        if (present(x2)) vect(vectInd(2)) = x2
        if (present(x3)) vect(vectInd(3)) = x3
        if (present(x4)) vect(vectInd(4)) = x4
        if (present(x5)) vect(vectInd(5)) = x5
        if (present(x6)) vect(vectInd(6)) = x6
        if (present(x7)) vect(vectInd(7)) = x7
        if (present(x8)) vect(vectInd(8)) = x8
        if (present(x9)) vect(vectInd(9)) = x9
        if (present(x10)) vect(vectInd(10)) = x10
        if (present(x11)) vect(vectInd(11)) = x11
        if (present(x12)) vect(vectInd(12)) = x12
        if (present(x13)) vect(vectInd(13)) = x13
        if (present(x14)) vect(vectInd(14)) = x14
        if (present(x15)) vect(vectInd(15)) = x15
        if (present(x16)) vect(vectInd(16)) = x16
        if (present(x17)) vect(vectInd(17)) = x17
    else
        ASSERT(.false.)
    end if

end subroutine
