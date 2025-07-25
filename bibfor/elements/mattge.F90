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
subroutine mattge(nomte, dtild, sina, cosa, r, &
                  jacp, vf, dfds, rtangi)
    implicit none
!
#include "asterfort/btkb.h"
#include "blas/dscal.h"
    character(len=16) :: nomte
    real(kind=8) :: sina, cosa, r, jacp, vf(*), dfds(*)
    real(kind=8) :: mats(5, 9), mats1(3, 9), dtild(5, 5), dtild1(3, 3)
    real(kind=8) :: rtangi(9, 9)
    real(kind=8) :: dtilds(5, 9), dtildt(3, 9)
!
!
!     CALCULS DE LA MATRICE TANGENTE
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j
    blas_int :: b_incx, b_n
!-----------------------------------------------------------------------
    if (nomte .eq. 'MECXSE3') then
!
        mats(1, 1) = -sina*dfds(1)
        mats(1, 2) = cosa*dfds(1)
        mats(1, 3) = 0.d0
        mats(1, 4) = -sina*dfds(2)
        mats(1, 5) = cosa*dfds(2)
        mats(1, 6) = 0.d0
        mats(1, 7) = -sina*dfds(3)
        mats(1, 8) = cosa*dfds(3)
        mats(1, 9) = 0.d0
!
        mats(2, 1) = 0.d0
        mats(2, 2) = 0.d0
        mats(2, 3) = dfds(1)
        mats(2, 4) = 0.d0
        mats(2, 5) = 0.d0
        mats(2, 6) = dfds(2)
        mats(2, 7) = 0.d0
        mats(2, 8) = 0.d0
        mats(2, 9) = dfds(3)
!
        mats(3, 1) = vf(1)/r
        mats(3, 2) = 0.d0
        mats(3, 3) = 0.d0
        mats(3, 4) = vf(2)/r
        mats(3, 5) = 0.d0
        mats(3, 6) = 0.d0
        mats(3, 7) = vf(3)/r
        mats(3, 8) = 0.d0
        mats(3, 9) = 0.d0
!
        mats(4, 1) = 0.d0
        mats(4, 2) = 0.d0
        mats(4, 3) = -sina*vf(1)/r
        mats(4, 4) = 0.d0
        mats(4, 5) = 0.d0
        mats(4, 6) = -sina*vf(2)/r
        mats(4, 7) = 0.d0
        mats(4, 8) = 0.d0
        mats(4, 9) = -sina*vf(3)/r
!
        mats(5, 1) = cosa*dfds(1)
        mats(5, 2) = sina*dfds(1)
        mats(5, 3) = vf(1)
        mats(5, 4) = cosa*dfds(2)
        mats(5, 5) = sina*dfds(2)
        mats(5, 6) = vf(2)
        mats(5, 7) = cosa*dfds(3)
        mats(5, 8) = sina*dfds(3)
        mats(5, 9) = vf(3)
!
        b_n = to_blas_int(25)
        b_incx = to_blas_int(1)
        call dscal(b_n, jacp, dtild, b_incx)
!
        call btkb(5, 9, 9, dtild, mats, &
                  dtilds, rtangi)
!
    else
!
        do i = 1, 3
            do j = 1, 3
                dtild1(i, j) = dtild(i, j)
            end do
        end do
!
        mats1(1, 1) = -sina*dfds(1)
        mats1(1, 2) = cosa*dfds(1)
        mats1(1, 3) = 0.d0
        mats1(1, 4) = -sina*dfds(2)
        mats1(1, 5) = cosa*dfds(2)
        mats1(1, 6) = 0.d0
        mats1(1, 7) = -sina*dfds(3)
        mats1(1, 8) = cosa*dfds(3)
        mats1(1, 9) = 0.d0
!
        mats1(2, 1) = 0.d0
        mats1(2, 2) = 0.d0
        mats1(2, 3) = dfds(1)
        mats1(2, 4) = 0.d0
        mats1(2, 5) = 0.d0
        mats1(2, 6) = dfds(2)
        mats1(2, 7) = 0.d0
        mats1(2, 8) = 0.d0
        mats1(2, 9) = dfds(3)
!
        mats1(3, 1) = cosa*dfds(1)
        mats1(3, 2) = sina*dfds(1)
        mats1(3, 3) = vf(1)
        mats1(3, 4) = cosa*dfds(2)
        mats1(3, 5) = sina*dfds(2)
        mats1(3, 6) = vf(2)
        mats1(3, 7) = cosa*dfds(3)
        mats1(3, 8) = sina*dfds(3)
        mats1(3, 9) = vf(3)
!
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, jacp, dtild1, b_incx)
!
        call btkb(3, 9, 9, dtild1, mats1, &
                  dtildt, rtangi)
!
    end if
!
end subroutine
