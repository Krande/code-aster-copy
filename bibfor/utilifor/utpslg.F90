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

subroutine utpslg(nn, nc, p, sl, sg)
    implicit none
!
#include "asterfort/mavec.h"
#include "asterfort/vecma.h"
!
    integer(kind=8)      :: nn, nc
    real(kind=8) :: p(3, 3), sl(*), sg(*)
!
! --------------------------------------------------------------------------------------------------
!
!     PASSAGE EN 3D D'UNE MATRICE TRIANGULAIRE DE NN*NC LIGNES
!     DU REPERE LOCAL AU REPERE GLOBAL
!
! --------------------------------------------------------------------------------------------------
!
!   In
!       nn   nombre de noeuds
!       nc   nombre de composantes
!       p    matrice de passage 3d de global a local
!       sl   nn*nc composantes de la triangulaire sl dans local
!   Out
!       sg   nn*nc composantes de la triangulaire sg dans global
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)      :: in(3)
!
    real(kind=8) :: r(9)
    real(kind=8) :: ml14(14, 14), mr14(14, 14), mtr14(14, 14), mv14(14, 14)
    real(kind=8) :: ml16(16, 16), mr16(16, 16), mtr16(16, 16), mv16(16, 16)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, l, m, n, nb
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
                    ! bloc diagonal
                    r(1) = sl(in(1)+1)
                    r(2) = sl(in(2)+1)
                    r(3) = sl(in(3)+1)
                    r(4) = sl(in(2)+1)
                    r(5) = sl(in(2)+2)
                    r(6) = sl(in(3)+2)
                    r(7) = sl(in(3)+1)
                    r(8) = sl(in(3)+2)
                    r(9) = sl(in(3)+3)
                    do m = 1, 3
                        do n = 1, m
                            sg(in(m)+n) = 0.0
                            do l = 1, 3
                                sg(in(m)+n) = sg(in(m)+n)+p(l, m)*(r(3*(l-1)+1)*p(1, n)+ &
                                                                   r(3*(l-1)+2)*p(2, n)+ &
                                                                   r(3*(l-1)+3)*p(3, n))
                            end do
                        end do
                    end do
                else
                    ! bloc extra - diagonal
                    do m = 1, 3
                        do n = 1, 3
                            sg(in(m)+n) = 0.0
                            do l = 1, 3
                                sg(in(m)+n) = sg(in(m)+n)+p(l, m)*(sl(in(l)+1)*p(1, n)+ &
                                                                   sl(in(l)+2)*p(2, n)+ &
                                                                   sl(in(l)+3)*p(3, n))
                            end do
                        end do
                    end do
                end if
            end do
        end do
!
    else if (mod(nc, 3) .eq. 1) then
        mr14(:, :) = 0.0
        do i = 1, 3
            do j = 1, 3
                mr14(i, j) = p(i, j)
                mr14(i+3, j+3) = p(i, j)
                mr14(i+7, j+7) = p(i, j)
                mr14(i+10, j+10) = p(i, j)
            end do
        end do
        mr14(7, 7) = 1.0
        mr14(14, 14) = 1.0
        mtr14 = transpose(mr14)
        call vecma(sl, 105, ml14, 14)
        mv14 = matmul(mtr14, ml14)
        mtr14 = matmul(mv14, mr14)
        call mavec(mtr14, 14, sg, 105)
!
    else if (mod(nc, 3) .eq. 2) then
        mr16(:, :) = 0.0
        do i = 1, 3
            do j = 1, 3
                mr16(i, j) = p(i, j)
                mr16(i+3, j+3) = p(i, j)
                mr16(i+8, j+8) = p(i, j)
                mr16(i+11, j+11) = p(i, j)
            end do
        end do
        mr16(7, 7) = 1.0
        mr16(8, 8) = 1.0
        mr16(15, 15) = 1.0
        mr16(16, 16) = 1.0
        mtr16 = transpose(mr16)
        call vecma(sl, 136, ml16, 16)
        mv16 = matmul(mtr16, ml16)
        mtr16 = matmul(mv16, mr16)
        call mavec(mtr16, 16, sg, 136)
    end if
!
end subroutine
