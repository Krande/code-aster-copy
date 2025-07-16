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

subroutine ut2mgl(nn, nc, p, sg, sl)
    implicit none
#include "asterfort/mavec.h"
#include "asterfort/vecma.h"
    real(kind=8) :: p(3, 3), sg(*), sl(*)
!     ------------------------------------------------------------------
!     PASSAGE EN 2D D'UNE MATRICE TRIANGULAIRE DE NN*NC LIGNES
!     DU REPERE GLOBAL AU REPERE LOCAL
!     ------------------------------------------------------------------
!IN   I   NN   NOMBRE DE NOEUDS
!IN   I   NC   NOMBRE DE COMPOSANTES
!IN   R   P    MATRICE DE PASSAGE 3D DE GLOBAL A LOCAL
!IN   R   SG   NN*NC COMPOSANTES DE LA TRIANGULAIRE SG DANS GLOBAL
!OUT  R   SL   NN*NC COMPOSANTES DE LA TRIANGULAIRE SL DANS LOCALL
!     ------------------------------------------------------------------
    real(kind=8) :: r(4)
    real(kind=8) :: ml6(6, 6), mr6(6, 6), mtr6(6, 6), mv6(6, 6)
    real(kind=8) :: ml3(3, 3), mr3(3, 3), mtr3(3, 3), mv3(3, 3)
    integer(kind=8) :: in(2)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, l, m, n, nb
    integer(kind=8) :: nc, nn
    real(kind=8) :: un, zero
!-----------------------------------------------------------------------
    zero = 0.0d0
    un = 1.0d0
!
    if (mod(nc, 2) .eq. 0) then
        nb = nn*nc/2
        do i = 1, nb
            k = 2*(i-1)
            do j = 1, i
                in(1) = k*(k+1)/2+2*(j-1)
                in(2) = (k+1)*(k+2)/2+2*(j-1)
                if (i .eq. j) then
!             --------- BLOC DIAGONAL
                    r(1) = sg(in(1)+1)
                    r(2) = sg(in(2)+1)
                    r(3) = sg(in(2)+1)
                    r(4) = sg(in(2)+2)
                    do m = 1, 2
                        do n = 1, m
                            sl(in(m)+n) = zero
                            do l = 1, 2
                                sl(in(m)+n) = sl(in(m)+n)+p(m, l)* &
                                              (r(2*(l-1)+1)*p(n, 1)+r(2*(l-1)+2)*p(n, 2))
                            end do
                        end do
                    end do
                else
!             --------- BLOC EXTRA - DIAGONAL
                    do m = 1, 2
                        do n = 1, 2
                            sl(in(m)+n) = zero
                            do l = 1, 2
                                sl(in(m)+n) = sl(in(m)+n)+p(m, l)* &
                                              (sg(in(l)+1)*p(n, 1)+sg(in(l)+2)*p(n, 2))
                            end do
                        end do
                    end do
                end if
            end do
        end do

!
    else if (mod(nc, 2) .eq. 1) then
        if (nc*nn .eq. 3) then
            mtr3(:, :) = zero
!
            mtr3(1, 1) = p(1, 1)
            mtr3(1, 2) = p(1, 2)
            mtr3(2, 1) = p(2, 1)
            mtr3(2, 2) = p(2, 2)
            mtr3(3, 3) = un
!
            mr3 = transpose(mtr3)
            call vecma(sg, 6, ml3, 3)
            mv3 = matmul(mtr3, ml3)
            mtr3 = matmul(mv3, mr3)
            call mavec(mtr3, 3, sl, 6)
!
        else if (nc*nn .eq. 6) then
            mtr6(:, :) = zero
!
            mtr6(1, 1) = p(1, 1)
            mtr6(1, 2) = p(1, 2)
            mtr6(2, 1) = p(2, 1)
            mtr6(2, 2) = p(2, 2)
            mtr6(3, 3) = un
!
            mtr6(4, 4) = p(1, 1)
            mtr6(4, 5) = p(1, 2)
            mtr6(5, 4) = p(2, 1)
            mtr6(5, 5) = p(2, 2)
            mtr6(6, 6) = un
!
            mr6 = transpose(mtr6)
            call vecma(sg, 21, ml6, 6)
            mv6 = matmul(mtr6, ml6)
            mtr6 = matmul(mv6, mr6)
            call mavec(mtr6, 6, sl, 21)
!
        end if
!
    end if
!
end subroutine ut2mgl
