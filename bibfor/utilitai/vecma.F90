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
subroutine vecma(mv, n, mp, m)
    implicit none
#include "asterfort/assert.h"
    integer(kind=8) :: n, m
    real(kind=8) :: mv(n), mp(m, m)
!       PASSAGE DEMI-MATRICE COLONNE VECTEUR (N)  > MATRICE PLEINE (M*M)
!       IN      MV = VECTEUR DEMI - MATRICE STOCKE COLONNE , LONGUEUR N
!       OUT     MP = MATRICE PLEINE (M*M)
!       ----------------------------------------------------------------
    integer(kind=8) :: k, i, j
! ----------------------------------------------------------------------
!
    ASSERT(n .ge. m*(m+1)/2)
!
    k = 0
    do i = 1, m
        do j = 1, m
            k = k+1
            mp(i, j) = mv(k)
            mp(j, i) = mv(k)
            if (j .eq. i) goto 10
        end do
10      continue
    end do
!
end subroutine
