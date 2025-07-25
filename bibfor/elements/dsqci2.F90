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
subroutine dsqci2(qsi, eta, caraq4, hft2, hmft2, &
                  bcb, bca, bcm)
    implicit none
    real(kind=8) :: qsi, eta, caraq4(*), hft2(2, 6), hmft2(2, 6)
    real(kind=8) :: bcb(2, 12), bcm(2, 8), bca(2, 4)
!     MATRICES BCB(2,12), BCA(2,4), BCM(2,8) AU POINT QSI, ETA POUR DSQ
!     -----------------------------------------------------------------
    integer(kind=8) :: i, j, k
    real(kind=8) :: peta, meta, pqsi, mqsi, c(4), s(4)
    real(kind=8) :: tb(6, 12), ta(6, 4), tc(6, 8)
!     ------------------------------------------------------------------
!
    c(1) = caraq4(13)
    c(2) = caraq4(14)
    c(3) = caraq4(15)
    c(4) = caraq4(16)
    s(1) = caraq4(17)
    s(2) = caraq4(18)
    s(3) = caraq4(19)
    s(4) = caraq4(20)
!
    peta = 1.d0+eta
    meta = 1.d0-eta
    pqsi = 1.d0+qsi
    mqsi = 1.d0-qsi
!
    do k = 1, 6
        do j = 1, 12
            tb(k, j) = 0.d0
        end do
    end do
    tb(3, 2) = 0.25d0
    tb(3, 5) = -0.25d0
    tb(3, 8) = 0.25d0
    tb(3, 11) = -0.25d0
    tb(6, 3) = 0.25d0
    tb(6, 6) = -0.25d0
    tb(6, 9) = 0.25d0
    tb(6, 12) = -0.25d0
!
    do i = 1, 6
        do j = 1, 4
            ta(i, j) = 0.d0
        end do
    end do
    ta(1, 1) = -meta*c(1)
    ta(1, 3) = -peta*c(3)
    ta(2, 2) = -pqsi*c(2)
    ta(2, 4) = -mqsi*c(4)
    ta(3, 1) = qsi*c(1)
    ta(3, 2) = -eta*c(2)
    ta(3, 3) = -qsi*c(3)
    ta(3, 4) = eta*c(4)
    ta(4, 1) = -meta*s(1)
    ta(4, 3) = -peta*s(3)
    ta(5, 2) = -pqsi*s(2)
    ta(5, 4) = -mqsi*s(4)
    ta(6, 1) = qsi*s(1)
    ta(6, 2) = -eta*s(2)
    ta(6, 3) = -qsi*s(3)
    ta(6, 4) = eta*s(4)
!
    do i = 1, 6
        do j = 1, 8
            tc(i, j) = 0.d0
        end do
    end do
    tc(3, 1) = 0.25d0
    tc(6, 2) = 0.25d0
    tc(3, 3) = -0.25d0
    tc(6, 4) = -0.25d0
    tc(3, 5) = 0.25d0
    tc(6, 6) = 0.25d0
    tc(3, 7) = -0.25d0
    tc(6, 8) = -0.25d0
!
!     -------------- BCB = HFT2.TB -----------------------------------
    do i = 1, 2
        do j = 1, 12
            bcb(i, j) = 0.d0
        end do
    end do
    do j = 1, 12
        do k = 1, 6
            bcb(1, j) = bcb(1, j)+hft2(1, k)*tb(k, j)
            bcb(2, j) = bcb(2, j)+hft2(2, k)*tb(k, j)
        end do
    end do
!     -------------- BCA = HFT2.TA -----------------------------------
    do i = 1, 2
        do j = 1, 4
            bca(i, j) = 0.d0
        end do
    end do
    do j = 1, 4
        do k = 1, 6
            bca(1, j) = bca(1, j)+hft2(1, k)*ta(k, j)
            bca(2, j) = bca(2, j)+hft2(2, k)*ta(k, j)
        end do
    end do
!     -------------- BCM = HMFT2.TC ----------------------------------
    do i = 1, 2
        do j = 1, 8
            bcm(i, j) = 0.d0
        end do
    end do
    do j = 1, 8
        do k = 1, 6
            bcm(1, j) = bcm(1, j)+hmft2(1, k)*tc(k, j)
            bcm(2, j) = bcm(2, j)+hmft2(2, k)*tc(k, j)
        end do
    end do
!
end subroutine
