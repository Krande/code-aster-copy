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
subroutine pacou4(a, n, c, d, sing)
    implicit none
!
! ARGUMENTS
! ---------
#include "asterf_types.h"
#include "jeveux.h"
    integer(kind=8) :: n
    real(kind=8) :: a(n, *), c(*), d(*)
    aster_logical :: sing
!
! --------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    real(kind=8) :: scale, sigma, sum, tau
!-----------------------------------------------------------------------
    sing = .false.
    scale = 0.0d0
!
    do k = 1, n-1
!
        do i = k, n
            scale = max(scale, abs(a(i, k)))
        end do
!
        if (abs(scale) .le. 1.0d-30) then
            sing = .true.
            c(k) = 0.0d0
            d(k) = 0.0d0
!
        else
            do i = k, n
                a(i, k) = a(i, k)/scale
            end do
!
            sum = 0.0d0
            do i = k, n
                sum = sum+a(i, k)**2
            end do
!
            sigma = sign(sqrt(sum), a(k, k))
            a(k, k) = a(k, k)+sigma
            c(k) = sigma*a(k, k)
            d(k) = -scale*sigma
!
            do j = k+1, n
!
                sum = 0.0d0
                do i = k, n
                    sum = sum+a(i, k)*a(i, j)
                end do
!
                tau = sum/c(k)
                do i = k, n
                    a(i, j) = a(i, j)-tau*a(i, k)
                end do
            end do
!
        end if
    end do
!
    d(n) = a(n, n)
    if (abs(d(n)) .le. 1.0d-30) sing = .true.
!
end subroutine
