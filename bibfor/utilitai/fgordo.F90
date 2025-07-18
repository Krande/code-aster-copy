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
subroutine fgordo(nbextr, ext, ord)
    implicit none
#include "jeveux.h"
    real(kind=8) :: ext(*), ord(*)
    integer(kind=8) :: nbextr
!     RANGE LES EXTREMAS PAR AMPLITUDE DECROISSANTE
!     -----------------------------------------------------------------
! IN  NBEXTR : I   : NOMBRE  D'EXTREMUM DE LA FONCTION
! IN  EXT    : R   : VALEURS DES EXTREMA
! OUT ORD    : R   : VALEURS DES EXTREMA REORDONNES
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
!-----------------------------------------------------------------------
    if (ext(1) .lt. ext(2)) then
        ord(1) = ext(1)
        ord(2) = ext(2)
    else
        ord(1) = ext(2)
        ord(2) = ext(1)
    end if
!
    do i = 3, nbextr
        if (ext(i) .lt. ord(1)) then
            do k = i, 2, -1
                ord(k) = ord(k-1)
            end do
            ord(1) = ext(i)
            goto 1
        end if
        if (ext(i) .ge. ord(i-1)) then
            ord(i) = ext(i)
            goto 1
        end if
        do j = 1, i-2
            if ((ord(j) .le. ext(i)) .and. (ext(i) .lt. ord(j+1))) then
                do k = i, j+2, -1
                    ord(k) = ord(k-1)
                end do
                ord(j+1) = ext(i)
                goto 1
            end if
        end do
1       continue
    end do
!
end subroutine
