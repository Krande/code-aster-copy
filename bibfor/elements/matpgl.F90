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
subroutine matpgl(nb1, xr, plg)
    implicit none
!
    integer(kind=8) :: nb1
    real(kind=8) :: plg(9, 3, 3), xr(*)
!
!     CONSTRUCTION DE LA MATRICE DE PASSAGE GLOBAL LOCAL  PGL(9,3,3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ib, j, k, l
!-----------------------------------------------------------------------
    do ib = 1, nb1
        l = 9*(ib-1)
        do j = 1, 3
            do i = 1, 3
                k = l+(j-1)*3+i
                plg(ib, i, j) = xr(1090+k)
            end do
        end do
!
    end do
!
end subroutine
