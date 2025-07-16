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
subroutine ordis(listis, nbterm)
    implicit none
    integer(kind=8) :: listis(*), nbterm
!
!     REARRANGEMENT D'UNE LISTE D'ENTIERS PAR ORDRE CROISSANT
! ------------------------------------------------------------------
!
    integer(kind=8) :: i, j, k
!
    do j = 2, nbterm
        k = listis(j)
        do i = j-1, 1, -1
            if (listis(i) .le. k) goto 30
            listis(i+1) = listis(i)
        end do
        i = 0
30      continue
        listis(i+1) = k
    end do
!
end subroutine
