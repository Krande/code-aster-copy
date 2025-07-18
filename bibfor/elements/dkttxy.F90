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
subroutine dkttxy(codi, lcot, hft2, depf, vt)
    implicit none
    real(kind=8) :: codi(*), lcot(*), hft2(2, 6), depf(9), vt(2)
!     EFFORT TRANCHANT L'ELEMENT DE PLAQUE DKT
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, k, nno, nc
    real(kind=8) :: tkt(6, 9), cl(3), sl(3), cs(3), cu(3), su(3), bc(2, 9)
!     ------------------------------------------------------------------
!
    nno = 3
    nc = 3
!
    do i = 1, nc
        cl(i) = 6.d0*codi(i)/lcot(i)
        sl(i) = 6.d0*codi(nc+i)/lcot(i)
        cs(i) = 3.d0*codi(i)*codi(nc+i)
        cu(i) = 3.d0*codi(i)*codi(i)
        su(i) = 3.d0*codi(nc+i)*codi(nc+i)
    end do
!
    tkt(1, 1) = -cl(1)-cl(1)
    tkt(1, 2) = cu(1)+cu(1)
    tkt(1, 3) = cs(1)+cs(1)
    tkt(1, 4) = cl(1)+cl(1)
    tkt(1, 5) = cu(1)+cu(1)
    tkt(1, 6) = cs(1)+cs(1)
    tkt(1, 7) = 0.d0
    tkt(1, 8) = 0.d0
    tkt(1, 9) = 0.d0
    tkt(2, 1) = cl(3)+cl(3)
    tkt(2, 2) = cu(3)+cu(3)
    tkt(2, 3) = cs(3)+cs(3)
    tkt(2, 4) = 0.d0
    tkt(2, 5) = 0.d0
    tkt(2, 6) = 0.d0
    tkt(2, 7) = -cl(3)-cl(3)
    tkt(2, 8) = cu(3)+cu(3)
    tkt(2, 9) = cs(3)+cs(3)
    tkt(3, 1) = cl(3)-cl(1)
    tkt(3, 2) = cu(3)+cu(1)
    tkt(3, 3) = cs(3)+cs(1)
    tkt(3, 4) = cl(1)+cl(2)
    tkt(3, 5) = cu(1)-cu(2)
    tkt(3, 6) = cs(1)-cs(2)
    tkt(3, 7) = -cl(3)-cl(2)
    tkt(3, 8) = cu(3)-cu(2)
    tkt(3, 9) = cs(3)-cs(2)
    tkt(4, 1) = -sl(1)-sl(1)
    tkt(4, 2) = cs(1)+cs(1)
    tkt(4, 3) = su(1)+su(1)
    tkt(4, 4) = sl(1)+sl(1)
    tkt(4, 5) = cs(1)+cs(1)
    tkt(4, 6) = su(1)+su(1)
    tkt(4, 7) = 0.d0
    tkt(4, 8) = 0.d0
    tkt(4, 9) = 0.d0
    tkt(5, 1) = sl(3)+sl(3)
    tkt(5, 2) = cs(3)+cs(3)
    tkt(5, 3) = su(3)+su(3)
    tkt(5, 4) = 0.d0
    tkt(5, 5) = 0.d0
    tkt(5, 6) = 0.d0
    tkt(5, 7) = -sl(3)-sl(3)
    tkt(5, 8) = cs(3)+cs(3)
    tkt(5, 9) = su(3)+su(3)
    tkt(6, 1) = sl(3)-sl(1)
    tkt(6, 2) = cs(3)+cs(1)
    tkt(6, 3) = su(3)+su(1)
    tkt(6, 4) = sl(1)+sl(2)
    tkt(6, 5) = cs(1)-cs(2)
    tkt(6, 6) = su(1)-su(2)
    tkt(6, 7) = -sl(3)-sl(2)
    tkt(6, 8) = cs(3)-cs(2)
    tkt(6, 9) = su(3)-su(2)
!
!     ------ VT = HFT2.TKT.DEPF ---------------------------------------
!
    bc(:, :) = 0.d0
!
    do i = 1, 2
        do j = 1, 3*nno
            do k = 1, 6
                bc(i, j) = bc(i, j)+hft2(i, k)*tkt(k, j)
            end do
        end do
    end do
    vt(1) = 0.d0
    vt(2) = 0.d0
    do i = 1, 2
        do j = 1, 3*nno
            vt(i) = vt(i)+bc(i, j)*depf(j)
        end do
    end do
!
end subroutine
