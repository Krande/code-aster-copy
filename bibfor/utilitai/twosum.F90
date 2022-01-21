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

subroutine twosum(a, b, x, y)
  implicit none

#include "asterfort/assert.h"

  real(kind=8), intent(in) :: a, b
  real(kind=8), intent(out) :: x, y

! ------------------------------------------------------------------------------
!
! Perform the TwoSum algorithm x+y = a+b
! Error-free transformation of the sum of two floating point numbers
! ------------------------------------------------------------------------------
!
  real(kind=8) :: z

  x = a + b
  z = x - a
  y = ( a - ( x - z )) + ( b - z )

end subroutine twosum

