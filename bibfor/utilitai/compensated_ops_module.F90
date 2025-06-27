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

module compensated_ops_module

    implicit none
#include "asterc/fma_double.h"

    interface dot_product
        module procedure dot2FMA
    end interface dot_product

    interface sum
        module procedure sum2s
    end interface sum

    interface matmul
        module procedure matmul2
    end interface matmul

! -------------------------------------------------------
!
! Compensated algorithms from Ogita-Oishi-Rump
!
! -------------------------------------------------------

contains

    subroutine twosum(a, b, x, y)
        implicit none
        real(kind=8), intent(in) :: a, b
        real(kind=8), intent(out) :: x, y
        real(kind=8) :: z

        x = a+b
        z = x-a
        y = (a-(x-z))+(b-z)

    end subroutine twosum

    subroutine split(a, x, y)
        implicit none
        real(kind=8), intent(in) :: a
        real(kind=8), intent(out) :: x, y
        real(kind=8), parameter :: factor = 2**27+1
        real(kind=8) :: c

        c = factor*a
        x = c-(c-a)
        y = a-x

    end subroutine split

    subroutine twoproduct(a, b, x, y)
        implicit none
        real(kind=8), intent(in) :: a, b
        real(kind=8), intent(out) :: x, y
        real(kind=8) :: a1, a2, b1, b2

        x = a*b
        call split(a, a1, a2)
        call split(b, b1, b2)

        y = a2*b2-(((x-a1*b1)-a2*b1)-a1*b2)

    end subroutine twoproduct

    subroutine twoproductFMA(a, b, x, y)
        implicit none
        real(kind=8), intent(in) :: a, b
        real(kind=8), intent(out) :: x, y

        x = a*b
        y = fma_double(a, b, -x)

    end subroutine twoproductFMA

    function sum2s(vt) result(t)
        implicit none
        real(kind=8), intent(in) :: vt(:)
        real(kind=8) :: si, x, y, t
        integer(kind=8) :: i

        t = vt(1)
        si = 0.d0
        do i = 2, size(vt)
            call twosum(t, vt(i), x, y)
            t = x
            si = si+y
        end do
        t = t+si

    end function sum2s

    function dot2s(vx, vy) result(t)
        implicit none
        real(kind=8), intent(in) :: vx(:)
        real(kind=8), intent(in) :: vy(size(vx))
        real(kind=8) :: p, s, h, r, q, pp, t
        integer(kind=8) i

        call twoproduct(vx(1), vy(1), p, s)
        do i = 2, size(vx)
            call twoproduct(vx(i), vy(i), h, r)
            pp = p
            call twosum(pp, h, p, q)
            s = s+(q+r)
        end do
        t = p+s

    end function dot2s

    function dot2FMA(vx, vy) result(t)
        implicit none
        real(kind=8), intent(in) :: vx(:)
        real(kind=8), intent(in) :: vy(size(vx))
        real(kind=8) :: p, s, h, r, q, pp, t
        integer(kind=8) i

        call twoproductFMA(vx(1), vy(1), p, s)
        do i = 2, size(vx)
            call twoproductFMA(vx(i), vy(i), h, r)
            pp = p
            call twosum(pp, h, p, q)
            s = s+(q+r)
        end do
        t = p+s

    end function dot2FMA

    function matmul2(A, B) result(C)
        implicit none
        real(kind=8), dimension(:, :), intent(in) :: A, B
        real(kind=8), dimension(size(A, 1), size(B, 2)) :: C
        integer(kind=8) :: i, j

        do i = 1, size(A, 1)
            do j = 1, size(B, 2)
                C(i, j) = dot2FMA(A(i, :), B(:, j))
            end do
        end do

    end function matmul2

end module compensated_ops_module
