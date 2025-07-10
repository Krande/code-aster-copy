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

subroutine iniVect0(l, x1, x2, x3, x4, x5, &
                    x6, x7, x8, x9, x10, x11, x12, x13, x14, &
                    x15, x16, x17, x18, x19)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!   initialisation des arguments à 0.d0
!-----------------------------------------------------------------------
    implicit none
    integer(kind=8), intent(in) :: l
    real(kind=8), intent(out) :: x1(l)
    real(kind=8), optional, intent(out) :: x2(l), x3(l), x4(l), x5(l), x6(l), x7(l), x8(l)
    real(kind=8), optional, intent(out) :: x9(l), x10(l), x11(l), x12(l), x13(l), x14(l)
    real(kind=8), optional, intent(out) :: x15(l), x16(l), x17(l), x18(l), x19(l)
!-----------------------------------------------------------------------
!
    x1(:) = 0.d0

    if (present(x2)) x2(:) = 0.d0
    if (present(x3)) x3(:) = 0.d0
    if (present(x4)) x4(:) = 0.d0
    if (present(x5)) x5(:) = 0.d0
    if (present(x6)) x6(:) = 0.d0
    if (present(x7)) x7(:) = 0.d0
    if (present(x8)) x8(:) = 0.d0
    if (present(x9)) x9(:) = 0.d0
    if (present(x10)) x10(:) = 0.d0
    if (present(x11)) x11(:) = 0.d0
    if (present(x12)) x12(:) = 0.d0
    if (present(x13)) x13(:) = 0.d0
    if (present(x14)) x14(:) = 0.d0
    if (present(x15)) x15(:) = 0.d0
    if (present(x16)) x16(:) = 0.d0
    if (present(x17)) x17(:) = 0.d0
    if (present(x18)) x18(:) = 0.d0
    if (present(x19)) x19(:) = 0.d0

end subroutine
