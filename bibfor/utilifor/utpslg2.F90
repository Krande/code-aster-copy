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

subroutine utpslg2(nn, nc, p, sl, sg)
    implicit none
!
#include "asterfort/assert.h"
!
    integer(kind=8)      :: nn, nc
    real(kind=8) :: p(3, 3), sl(*), sg(*)
!
! --------------------------------------------------------------------------------------------------
!
!     PASSAGE EN 3D D'UNE MATRICE PLEINE DE NN*NC LIGNES
!     DU REPERE LOCAL AU REPERE GLOBAL
!
! --------------------------------------------------------------------------------------------------
!
!   In
!       nn   nombre de noeuds
!       nc   nombre de composantes
!       p    matrice de passage 3d de global a local
!       sl   nn*nc composantes de la sl pleine dans local
!   Out
!       sg   nn*nc composantes de la sg pleine dans global
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)      :: in(3)
!
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j, l, m, n, nb
!
    ASSERT(mod(nc, 3) .eq. 0)
    nb = nn*nc/3
    do i = 1, nb
        do j = 1, nb
            in(1) = 3*nb*3*(i-1)+3*(j-1)
            in(2) = 3*nb*(3*(i-1)+1)+3*(j-1)
            in(3) = 3*nb*(3*(i-1)+2)+3*(j-1)
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
        end do
    end do
!
end subroutine
