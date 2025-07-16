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

subroutine getR6Mat6(mat, x1, x2, x3, x4, x5, &
                     x6, x7, x8, x9, x10, x11, x12, x13, x14, &
                     x15, x16, x17)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!   récuperation des valeurs d'un tableau
!-----------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in) :: mat(6, *)
    real(kind=8), intent(out) :: x1(6)
    real(kind=8), optional, intent(out) :: x2(6), x3(6), x4(6), x5(6)
    real(kind=8), optional, intent(out) :: x6(6), x7(6), x8(6), x9(6)
    real(kind=8), optional, intent(out) :: x10(6), x11(6), x12(6), x13(6)
    real(kind=8), optional, intent(out) :: x14(6), x15(6), x16(6), x17(6)
!-----------------------------------------------------------------------
!
    x1(:) = mat(:, 1)
    if (present(x2)) x2(:) = mat(:, 2)
    if (present(x3)) x3(:) = mat(:, 3)
    if (present(x4)) x4(:) = mat(:, 4)
    if (present(x5)) x5(:) = mat(:, 5)
    if (present(x6)) x6(:) = mat(:, 6)
    if (present(x7)) x7(:) = mat(:, 7)
    if (present(x8)) x8(:) = mat(:, 8)
    if (present(x9)) x9(:) = mat(:, 9)
    if (present(x10)) x10(:) = mat(:, 10)
    if (present(x11)) x11(:) = mat(:, 11)
    if (present(x12)) x12(:) = mat(:, 12)
    if (present(x13)) x13(:) = mat(:, 13)
    if (present(x14)) x14(:) = mat(:, 14)
    if (present(x15)) x15(:) = mat(:, 15)
    if (present(x16)) x16(:) = mat(:, 16)
    if (present(x17)) x17(:) = mat(:, 17)

end subroutine
