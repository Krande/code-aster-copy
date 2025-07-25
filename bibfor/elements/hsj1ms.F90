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
subroutine hsj1ms(epais, vectg, vectt, hsfm, hss, &
                  hsj1m, hsj1s)
    implicit none
#include "asterfort/jacbm1.h"
    real(kind=8) :: epais
    real(kind=8) :: vectg(2, 3), vectt(3, 3), matj(3, 3), jm1(3, 3), detj
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
!
!     CONSTRUCTION DE J-1 AUX PTS D'INTEGRATION REDUITS.   J = JACOBIEN
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, j1, jb, k, k1
!-----------------------------------------------------------------------
    call jacbm1(epais, vectg, vectt, matj, jm1, &
                detj)
!
!     CONSTRUCTION DE HFM * S * JTILD-1  ET  HS * S * JTILD-1
!                                            AUX PTS DE GAUSS REDUITS
!
!     PARTITION EN BLOCS DE HSFM ET HSS
!
!     HFM * S = HSFM = ( HSFM11 , HSFM12 , HSFM13 )  OU   HSFM1J (3,3)
!
!     HS  * S = HSS =  ( HSS11 , HSS12 , HSS13 )     OU   HSS1J (2,3)
!
!               J-1 0  0
!     JTILD-1 = 0  J-1 0                        OU   J-1   (3,3)
!               0   0  J-1
!
!           IB=1,1
    do jb = 1, 3
        do j = 1, 3
            do i = 1, 2
                j1 = j+3*(jb-1)
                hsj1m(i, j1) = 0.d0
                hsj1s(i, j1) = 0.d0
                do k = 1, 3
                    k1 = k+3*(jb-1)
                    hsj1s(i, j1) = hsj1s(i, j1)+hss(i, k1)*jm1(k, j)
                    hsj1m(i, j1) = hsj1m(i, j1)+hsfm(i, k1)*jm1(k, j)
                end do
            end do
            j1 = j+3*(jb-1)
            hsj1m(3, j1) = 0.d0
            do k = 1, 3
                k1 = k+3*(jb-1)
                hsj1m(3, j1) = hsj1m(3, j1)+hsfm(3, k1)*jm1(k, j)
            end do
        end do
    end do
!
end subroutine
