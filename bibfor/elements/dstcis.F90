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
subroutine dstcis(dci, carat3, hft2, bca, an)
    implicit none
#include "asterfort/mgauss.h"
    real(kind=8) :: dci(2, 2), carat3(*), hft2(2, 6), bca(2, 3), an(3, 9)
!     MATRICES BCA(2,3) ET AN(3,9) DU CISAILLEMENT POUR LE DST
!     --------------------------------------------------------
!
    integer(kind=8) :: i, j, k, iret
    real(kind=8) :: l(3), c(3), s(3), x(3), y(3), det
    real(kind=8) :: ta(6, 3), db(2, 3), aa(3, 3), aai(3, 3), aw(3, 9)
!     ------------------------------------------------------------------
    c(1) = carat3(16)
    c(2) = carat3(17)
    c(3) = carat3(18)
    s(1) = carat3(19)
    s(2) = carat3(20)
    s(3) = carat3(21)
    l(1) = carat3(13)
    l(2) = carat3(14)
    l(3) = carat3(15)
    x(1) = carat3(1)
    x(2) = carat3(2)
    x(3) = carat3(3)
    y(1) = carat3(4)
    y(2) = carat3(5)
    y(3) = carat3(6)
!
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
!     -------------- PRODUIT HFT2.TA -----------------------------------
    do i = 1, 2
        do j = 1, 3
            bca(i, j) = 0.d0
        end do
    end do
    do j = 1, 3
        do k = 1, 6
            bca(1, j) = bca(1, j)+hft2(1, k)*ta(k, j)
            bca(2, j) = bca(2, j)+hft2(2, k)*ta(k, j)
        end do
    end do
!     -------------- PRODUIT DCI.BCA -----------------------------------
    do j = 1, 3
        db(1, j) = dci(1, 1)*bca(1, j)+dci(1, 2)*bca(2, j)
        db(2, j) = dci(2, 1)*bca(1, j)+dci(2, 2)*bca(2, j)
    end do
!     -------------- CALCUL DE AA --------------------------------------
    do i = 1, 3
        do j = 1, 3
            aa(i, j) = -(x(i)*db(1, j)+y(i)*db(2, j))
        end do
        aa(i, i) = aa(i, i)+2.d0/3.d0*l(i)
    end do
!     -------------- INVERSION DE AA -----------------------------------
    do i = 1, 3
        do j = 1, 3
            aai(i, j) = 0.d0
        end do
    end do
    do i = 1, 3
        aai(i, i) = 1.d0
    end do
    call mgauss('NFVP', aa, aai, 3, 3, &
                3, det, iret)
!
!     -------------- CALCUL DE AW --------------------------------------
    do i = 1, 3
        do j = 1, 9
            aw(i, j) = 0.d0
        end do
    end do
    aw(1, 1) = 1.d0
    aw(1, 2) = -x(1)/2.d0
    aw(1, 3) = -y(1)/2.d0
    aw(1, 4) = -1.d0
    aw(1, 5) = -x(1)/2.d0
    aw(1, 6) = -y(1)/2.d0
    aw(2, 4) = 1.d0
    aw(2, 5) = -x(2)/2.d0
    aw(2, 6) = -y(2)/2.d0
    aw(2, 7) = -1.d0
    aw(2, 8) = -x(2)/2.d0
    aw(2, 9) = -y(2)/2.d0
    aw(3, 1) = -1.d0
    aw(3, 2) = -x(3)/2.d0
    aw(3, 3) = -y(3)/2.d0
    aw(3, 7) = 1.d0
    aw(3, 8) = -x(3)/2.d0
    aw(3, 9) = -y(3)/2.d0
!
!     -------------- PRODUIT AAI.AW ------------------------------------
    do i = 1, 3
        do j = 1, 9
            an(i, j) = 0.d0
        end do
    end do
    do i = 1, 3
        do k = 1, 3
            do j = 1, 9
                an(i, j) = an(i, j)+aai(i, k)*aw(k, j)
            end do
        end do
    end do
!
end subroutine
