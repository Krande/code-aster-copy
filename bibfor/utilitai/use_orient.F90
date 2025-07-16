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

! aslint: disable=W0413
function use_orient(vect, size, epsilon)
    implicit none
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
!
!   Tell if the vector defines valid nautic angles.
!   The values must not be nul (less than 'epsilon'), Nan or undef.
!
!   Arguments:
!       vect: vector of nautic angles values.
!       size: size of the vector.
!       epsilon: precision, default to r8prem()
!   Returns:
!       True if the angles are really defined.
!
    aster_logical :: use_orient
    integer(kind=8), intent(in) :: size
    real(kind=8), intent(in) :: vect(size)
    real(kind=8), intent(in), optional :: epsilon
!
    aster_logical :: bool
    integer(kind=8) :: i
    real(kind=8) :: eps

    if (present(epsilon)) then
        eps = epsilon
    else
        eps = r8prem()
    end if

    bool = ASTER_FALSE
    do i = 1, size
        if (isnan(vect(i)) .or. vect(i) .eq. r8maem()) then
            exit
        end if
        if (abs(vect(i)) .gt. eps) then
            bool = ASTER_TRUE
            exit
        end if
    end do
    use_orient = bool

end function use_orient
