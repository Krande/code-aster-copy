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
subroutine dstlxy(codi, hlt2, an, depf, lambda)
    implicit none
    real(kind=8) :: hlt2(4, 6), an(3, 9), depf(9), codi(*), lambda(4)
!     'LAMBDA' DE L'ELEMENT DE PLAQUE DST
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    real(kind=8) :: c(3), s(3)
    real(kind=8) :: ta(6, 3), bla(4, 3), bln(4, 9)
!     ------------------------------------------------------------------
!
!       -------------- BLA = HLT2.TA ------------------------------
    c(1) = codi(1)
    c(2) = codi(2)
    c(3) = codi(3)
    s(1) = codi(4)
    s(2) = codi(5)
    s(3) = codi(6)
    do k = 1, 6
        do j = 1, 3
            ta(k, j) = 0.d0
        end do
    end do
    ta(1, 1) = -8.d0*c(1)
    ta(2, 3) = -8.d0*c(3)
    ta(3, 1) = -4.d0*c(1)
    ta(3, 2) = 4.d0*c(2)
    ta(3, 3) = -4.d0*c(3)
    ta(4, 1) = -8.d0*s(1)
    ta(5, 3) = -8.d0*s(3)
    ta(6, 1) = -4.d0*s(1)
    ta(6, 2) = 4.d0*s(2)
    ta(6, 3) = -4.d0*s(3)
    do i = 1, 4
        do j = 1, 3
            bla(i, j) = 0.d0
            do k = 1, 6
                bla(i, j) = bla(i, j)+hlt2(i, k)*ta(k, j)
            end do
        end do
    end do
!       -------- LAMBDA = BLA.AN.DEPF ------------------------------
    do i = 1, 4
        lambda(i) = 0.d0
    end do
    do i = 1, 4
        do j = 1, 9
            bln(i, j) = 0.d0
            do k = 1, 3
                bln(i, j) = bln(i, j)+bla(i, k)*an(k, j)
            end do
            lambda(i) = lambda(i)+bln(i, j)*depf(j)
        end do
    end do
!
end subroutine
