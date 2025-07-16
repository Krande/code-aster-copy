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

subroutine getMat33Tab(tab, x1, x2, x3, x4, x5, &
                       x6, x7, x8, x9, x10, x11, x12, x13, x14, &
                       x15, x16, x17)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!   récuperation des valeurs d'un tableau
!-----------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in) :: tab(3, 3, *)
    real(kind=8), intent(out) :: x1(3, 3)
    real(kind=8), optional, intent(out) :: x2(3, 3), x3(3, 3), x4(3, 3), x5(3, 3)
    real(kind=8), optional, intent(out) :: x6(3, 3), x7(3, 3), x8(3, 3), x9(3, 3)
    real(kind=8), optional, intent(out) :: x10(3, 3), x11(3, 3), x12(3, 3), x13(3, 3)
    real(kind=8), optional, intent(out) :: x14(3, 3), x15(3, 3), x16(3, 3), x17(3, 3)
!-----------------------------------------------------------------------
!
    x1(:, :) = tab(:, :, 1)
    if (present(x2)) x2(:, :) = tab(:, :, 2)
    if (present(x3)) x3(:, :) = tab(:, :, 3)
    if (present(x4)) x4(:, :) = tab(:, :, 4)
    if (present(x5)) x5(:, :) = tab(:, :, 5)
    if (present(x6)) x6(:, :) = tab(:, :, 6)
    if (present(x7)) x7(:, :) = tab(:, :, 7)
    if (present(x8)) x8(:, :) = tab(:, :, 8)
    if (present(x9)) x9(:, :) = tab(:, :, 9)
    if (present(x10)) x10(:, :) = tab(:, :, 10)
    if (present(x11)) x11(:, :) = tab(:, :, 11)
    if (present(x12)) x12(:, :) = tab(:, :, 12)
    if (present(x13)) x13(:, :) = tab(:, :, 13)
    if (present(x14)) x14(:, :) = tab(:, :, 14)
    if (present(x15)) x15(:, :) = tab(:, :, 15)
    if (present(x16)) x16(:, :) = tab(:, :, 16)
    if (present(x17)) x17(:, :) = tab(:, :, 17)

end subroutine
