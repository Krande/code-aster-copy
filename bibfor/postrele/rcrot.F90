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

subroutine rcrot(nbabsc, phi, vale, sigm)
    implicit none
    integer(kind=8) :: nbabsc, i
    real(kind=8) :: phi, vale(4, nbabsc), sigm(4, nbabsc)
!     CALCULER LA TENSEUR DES CONTRAINTES ECRITE AVEC LA BASE LOCALE D'UN SEGMENT AXISYMETRIQUE
!
!
!     ------------------------------------------------------------------
    real(kind=8) :: rot(2, 2), valr(2, 2), sigr(2, 2)
! DEB ------------------------------------------------------------------
    rot(1, 1) = cos(phi)
    rot(1, 2) = -sin(phi)
    rot(2, 1) = sin(phi)
    rot(2, 2) = cos(phi)
    do i = 1, nbabsc
        valr(1, 1) = vale(1, i)
        valr(1, 2) = vale(4, i)
        valr(2, 1) = vale(4, i)
        valr(2, 2) = vale(2, i)
        sigr = matmul(matmul(transpose(rot), valr), rot)
        sigm(1, i) = sigr(1, 1)
        sigm(2, i) = sigr(2, 2)
        sigm(3, i) = vale(3, i)
        sigm(4, i) = sigr(1, 2)
    end do
!
!
end subroutine
