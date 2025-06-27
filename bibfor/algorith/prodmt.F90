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
subroutine prodmt(matsym, tens4, mat)
    implicit none
!     Product matric * tenseur4
! ----------------------------------------------------------------------
!   In matsym : matrix (XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ)
!   In tens4  : tensor of order 4
!   Out Mat = matsym * tens4
! ----------------------------------------------------------------------
!
    real(kind=8), intent(in) :: matsym(6), tens4(3, 3, 3, 3)
    real(kind=8), intent(out) :: mat(3, 3)
    real(kind=8) :: matS(3, 3)
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8) :: i, j, k, l
! ----------------------------------------------------------------------
!
    matS = 0.d0
    !
    matS(1, 1) = matsym(1)
    matS(2, 2) = matsym(2)
    matS(3, 3) = matsym(3)
    !
    matS(1, 2) = matsym(4)/rac2
    matS(2, 1) = matS(1, 2)
    matS(1, 3) = matsym(5)/rac2
    matS(3, 1) = matS(1, 3)
    matS(2, 3) = matsym(6)/rac2
    matS(3, 2) = matS(2, 3)
!
    mat = 0.d0
    do i = 1, 3
        do j = 1, 3
            do k = 1, 3
                do l = 1, 3
                    mat(i, j) = mat(i, j)+matS(k, l)*tens4(k, l, i, j)
                end do
            end do
        end do
    end do
!
end subroutine
