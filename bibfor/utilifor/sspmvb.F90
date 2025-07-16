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
subroutine sspmvb(n, m, mat, ad, t1, &
                  y)
!
!    VERSION FORTRAN DE LA ROUTINE CAL SSPMVA
!
! DEB ------------------------------------------------------------------
    implicit none
    integer(kind=8) :: n, m, ad(*)
    real(kind=8) :: mat(*), t1(*), y(*)
    integer(kind=8) :: i, j, jmin, jmax, rest, k0, k1, k2, k3, k4, k5, k6, k7
    rest = m-(m/8)*8
    if (rest .le. 3) then
!
        if (rest .le. 1) then
!
            if (rest .eq. 1) then
                k1 = ad(1)
                do i = 1, n
                    y(i) = y(i)-t1(1)*mat(k1)
                    k1 = k1+1
                end do
            end if
        else
            if (rest .eq. 2) then
                k1 = ad(1)
                k2 = ad(2)
                do i = 1, n
                    y(i) = y(i)-t1(1)*mat(k1)-t1(2)*mat(k2)
                    k1 = k1+1
                    k2 = k2+1
                end do
            else
                k1 = ad(1)
                k2 = ad(2)
                k3 = ad(3)
                do i = 1, n
                    y(i) = y(i)-t1(1)*mat(k1)-t1(2)*mat(k2)-t1(3)*mat(k3)
                    k1 = k1+1
                    k2 = k2+1
                    k3 = k3+1
                end do
            end if
        end if
    else
        if (rest .le. 5) then
!
            if (rest .eq. 4) then
                k1 = ad(1)
                k2 = ad(2)
                k3 = ad(3)
                k4 = ad(4)
                do i = 1, n
                    y(i) = y(i)-t1(1)*mat(k1)-t1(2)*mat(k2)-t1(3)*mat(k3)-t1(4)*mat(k4)
                    k1 = k1+1
                    k2 = k2+1
                    k3 = k3+1
                    k4 = k4+1
                end do
            else
                k1 = ad(1)
                k2 = ad(2)
                k3 = ad(3)
                k4 = ad(4)
                k5 = ad(5)
                do i = 1, n
                    y(i) = y(i)-t1(1)*mat(k1)-t1(2)*mat(k2)-t1(3)*mat(k3)-t1(4)*mat(k4)-&
                           & t1(5)*mat(k5)
                    k1 = k1+1
                    k2 = k2+1
                    k3 = k3+1
                    k4 = k4+1
                    k5 = k5+1
                end do
            end if
        else
            if (rest .eq. 6) then
                k1 = ad(1)
                k2 = ad(2)
                k3 = ad(3)
                k4 = ad(4)
                k5 = ad(5)
                k6 = ad(6)
                do i = 1, n
                    y(i) = y(i)-t1(1)*mat(k1)-t1(2)*mat(k2)-t1(3)*mat(k3)-t1(4)*mat(k4)-&
                           & t1(5)*mat(k5)-t1(6)*mat(k6)
                    k1 = k1+1
                    k2 = k2+1
                    k3 = k3+1
                    k4 = k4+1
                    k5 = k5+1
                    k6 = k6+1
                end do
!
            else
                k1 = ad(1)
                k2 = ad(2)
                k3 = ad(3)
                k4 = ad(4)
                k5 = ad(5)
                k6 = ad(6)
                k7 = ad(7)
                do i = 1, n
                    y(i) = y(i)-t1(1)*mat(k1)-t1(2)*mat(k2)-t1(3)*mat(k3)-t1(4)*mat(k4)-&
                           & t1(5)*mat(k5)-t1(6)*mat(k6)-t1(7)*mat(k7)
                    k1 = k1+1
                    k2 = k2+1
                    k3 = k3+1
                    k4 = k4+1
                    k5 = k5+1
                    k6 = k6+1
                    k7 = k7+1
                end do
            end if
        end if
    end if
    jmin = rest+8
    jmax = m
    if (jmax .ge. jmin) then
        do j = jmin, jmax, 8
            k0 = ad(j)
            k1 = ad(j-1)
            k2 = ad(j-2)
            k3 = ad(j-3)
            k4 = ad(j-4)
            k5 = ad(j-5)
            k6 = ad(j-6)
            k7 = ad(j-7)
            do i = 1, n
                y(i) = y(i)-t1(j)*mat(k0)-t1(j-1)*mat(k1)-t1(j-2)*mat(k2)-t1(j-3)*mat(k3&
                       &)-t1(j-4)*mat(k4)-t1(j-5)*mat(k5)-t1(j-6)*mat(k6)-t1(j-7)*mat(k7)
                k0 = k0+1
                k1 = k1+1
                k2 = k2+1
                k3 = k3+1
                k4 = k4+1
                k5 = k5+1
                k6 = k6+1
                k7 = k7+1
            end do
!
        end do
    end if
! FIN ------------------------------------------------------------------
end subroutine
