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
subroutine predia(a, b, diag, nno)
    implicit none
    integer(kind=8) :: i, ic, j, nno
!-----------------------------------------------------------------------
!
!    ESTIMATEUR ZZ (2-EME VERSION 92)
!
!      PRECONDITIONNEMENT PAR LA DIAGONALE DU SYSTEME LOCAL
!
    real(kind=8) :: a(9, 9), b(9, 4), diag(9)
    do i = 1, nno
        diag(i) = 1.d0/sqrt(a(i, i))
    end do
    do i = 1, nno
        do j = 1, i
            a(i, j) = a(i, j)*diag(i)*diag(j)
        end do
    end do
    do ic = 1, 4
        do i = 1, nno
            b(i, ic) = b(i, ic)*diag(i)
        end do
    end do
    do i = 1, nno
        do j = i+1, nno
            a(i, j) = a(j, i)
        end do
    end do
end subroutine
