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
subroutine utpsgl(nn, nc, p, sg, sl)
    implicit none
#include "asterfort/mavec.h"
#include "asterfort/vecma.h"
    real(kind=8) :: p(3, 3), sg(*), sl(*)
!     ------------------------------------------------------------------
!     PASSAGE EN 3D D'UNE MATRICE TRIANGULAIRE DE NN*NC LIGNES
!     DU REPERE GLOBAL AU REPERE LOCAL
!     ------------------------------------------------------------------
!IN   I   NN   NOMBRE DE NOEUDS
!IN   I   NC   NOMBRE DE COMPOSANTES
!IN   R   P    MATRICE DE PASSAGE 3D DE GLOBAL A LOCAL
!IN   R   SG   NN*NC COMPOSANTES DE LA TRIANGULAIRE SG DANS GLOBAL
!OUT  R   SL   NN*NC COMPOSANTES DE LA TRIANGULAIRE SL DANS LOCAL
!     ------------------------------------------------------------------
    real(kind=8) :: r(9), zero
    real(kind=8) :: ml14(14, 14), mr14(14, 14), mtr14(14, 14), mv14(14, 14)
    real(kind=8) :: ml16(16, 16), mr16(16, 16), mtr16(16, 16), mv16(16, 16)
    integer(kind=8) :: in(3)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, l, m, n, nb
    integer(kind=8) :: nc, nn
!-----------------------------------------------------------------------
    data zero/0.d0/
!
    if (mod(nc, 3) .eq. 0) then
        nb = nn*nc/3
        do i = 1, nb
            k = 3*(i-1)
            do j = 1, i
                in(1) = k*(k+1)/2+3*(j-1)
                in(2) = (k+1)*(k+2)/2+3*(j-1)
                in(3) = (k+2)*(k+3)/2+3*(j-1)
                if (i .eq. j) then
!             --------- BLOC DIAGONAL
                    r(1) = sg(in(1)+1)
                    r(2) = sg(in(2)+1)
                    r(3) = sg(in(3)+1)
                    r(4) = sg(in(2)+1)
                    r(5) = sg(in(2)+2)
                    r(6) = sg(in(3)+2)
                    r(7) = sg(in(3)+1)
                    r(8) = sg(in(3)+2)
                    r(9) = sg(in(3)+3)
                    do m = 1, 3
                        do n = 1, m
                            sl(in(m)+n) = zero
                            do l = 1, 3
                                sl(in(m)+n) = sl( &
                                              in(m)+n)+p(m, &
                                                         l)*( &
                                              r( &
                                              3*(l-1)+1)*p(n, 1)+r(3*(l-1)+2)*p(n, &
                                                                                2)+r(3*(l-1)+3 &
                                                                                     )*p(n, 3 &
                                                                                         ) &
                                              )
                            end do
                        end do
                    end do
                else
!             --------- BLOC EXTRA - DIAGONAL
                    do m = 1, 3
                        do n = 1, 3
                            sl(in(m)+n) = zero
                            do l = 1, 3
                                sl(in(m)+n) = sl( &
                                              in(m)+n)+p(m, &
                                                         l)*( &
                                              sg( &
                                              in(l)+1)*p(n, 1)+sg(in(l)+2)*p(n, &
                                                                             2)+sg(in(l)+3 &
                                                                                   )*p(n, 3 &
                                                                                       ) &
                                              )
                            end do
                        end do
                    end do
                end if
            end do
        end do
!
    else if (mod(nc, 3) .eq. 1) then
        do i = 1, 14
            do j = 1, 14
                mtr14(i, j) = 0.d0
            end do
        end do
        do i = 1, 3
            do j = 1, 3
                mtr14(i, j) = p(i, j)
                mtr14(i+3, j+3) = p(i, j)
                mtr14(i+7, j+7) = p(i, j)
                mtr14(i+10, j+10) = p(i, j)
            end do
        end do
        mtr14(7, 7) = 1.d0
        mtr14(14, 14) = 1.d0
        mr14 = transpose(mtr14)
        call vecma(sl, 105, ml14, 14)
        mv14 = matmul(mtr14, ml14)
        mtr14 = matmul(mv14, mr14)
        call mavec(mtr14, 14, sg, 105)
!
    else if (mod(nc, 3) .eq. 2) then
        do i = 1, 16
            do j = 1, 16
                mtr16(i, j) = 0.d0
            end do
        end do
        do i = 1, 3
            do j = 1, 3
                mtr16(i, j) = p(i, j)
                mtr16(i+3, j+3) = p(i, j)
                mtr16(i+8, j+8) = p(i, j)
                mtr16(i+11, j+11) = p(i, j)
            end do
        end do
        mtr16(7, 7) = 1.d0
        mtr16(8, 8) = 1.d0
        mtr16(15, 15) = 1.d0
        mtr16(16, 16) = 1.d0
        mr16 = transpose(mtr16)
        call vecma(sl, 136, ml16, 16)
        mv16 = matmul(mtr16, ml16)
        mtr16 = matmul(mv16, mr16)
        call mavec(mtr16, 16, sg, 136)
    end if
!
end subroutine
