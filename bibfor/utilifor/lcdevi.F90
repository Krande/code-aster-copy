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
subroutine lcdevi(a, d)
    implicit none
!       DEVIATEUR D UN TENSEUR (3X3) SOUS FORME VECTEUR  (6X1)
!       IN  A      :  TENSEUR
!       OUT D      :  DEVIATEUR DE A = A - 1/3 TR(A) I
!       ----------------------------------------------------------------
    integer(kind=8) :: ndt, ndi
    real(kind=8) :: a(6), d(6), ta
    common/tdim/ndt, ndi
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i
!-----------------------------------------------------------------------
    d = 0.d0
    ta = sum(a(1:ndi))/3.d0
!
    d(1:ndt) = a(1:ndt)
    d(1:ndi) = d(1:ndi)-ta

end subroutine
