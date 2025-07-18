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
subroutine pacou7(a, n, d, b)
    implicit none
!
! ARGUMENTS
! ---------
#include "jeveux.h"
    integer(kind=8) :: n
    real(kind=8) :: a(n, *), b(*), d(*)
! ---------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
    real(kind=8) :: sum
!-----------------------------------------------------------------------
    b(n) = b(n)/d(n)
!
    do i = n-1, 1, -1
!
        sum = 0.0d0
        do j = i+1, n
            sum = sum+a(i, j)*b(j)
        end do
!
        b(i) = (b(i)-sum)/d(i)
!
    end do
!
end subroutine
