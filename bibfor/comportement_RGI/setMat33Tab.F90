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

subroutine setMat33Tab(tab, x1, x2, x3, x4, x5, &
                       x6, x7, x8, x9, x10, x11, x12, x13, x14, &
                       x15, x16, x17)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!   récuperation des valeurs d'un tableau
!-----------------------------------------------------------------------
    implicit none
    real(kind=8), intent(inout) :: tab(3, 3, *)
    real(kind=8), intent(in) :: x1(3, 3)
    real(kind=8), optional, intent(in) :: x2(3, 3), x3(3, 3), x4(3, 3), x5(3, 3)
    real(kind=8), optional, intent(in) :: x6(3, 3), x7(3, 3), x8(3, 3), x9(3, 3)
    real(kind=8), optional, intent(in) :: x10(3, 3), x11(3, 3), x12(3, 3), x13(3, 3)
    real(kind=8), optional, intent(in) :: x14(3, 3), x15(3, 3), x16(3, 3), x17(3, 3)
!-----------------------------------------------------------------------
!
    tab(:, :, 1) = x1(:, :)
    if (present(x2)) tab(:, :, 2) = x2(:, :)
    if (present(x3)) tab(:, :, 3) = x3(:, :)
    if (present(x4)) tab(:, :, 4) = x4(:, :)
    if (present(x5)) tab(:, :, 5) = x5(:, :)
    if (present(x6)) tab(:, :, 6) = x6(:, :)
    if (present(x7)) tab(:, :, 7) = x7(:, :)
    if (present(x8)) tab(:, :, 8) = x8(:, :)
    if (present(x9)) tab(:, :, 9) = x9(:, :)
    if (present(x10)) tab(:, :, 10) = x10(:, :)
    if (present(x11)) tab(:, :, 11) = x11(:, :)
    if (present(x12)) tab(:, :, 12) = x12(:, :)
    if (present(x13)) tab(:, :, 13) = x13(:, :)
    if (present(x14)) tab(:, :, 14) = x14(:, :)
    if (present(x15)) tab(:, :, 15) = x15(:, :)
    if (present(x16)) tab(:, :, 16) = x16(:, :)
    if (present(x17)) tab(:, :, 17) = x17(:, :)

end subroutine
