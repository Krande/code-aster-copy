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

subroutine getIntVect(vect, x1, x2, x3, x4, x5, &
                      x6, x7, x8, x9, x10, x11, x12, x13, x14, &
                      x15, x16, x17)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!   récupèration de valeurs entières dans un tableau
!-----------------------------------------------------------------------
    implicit none
    integer(kind=8), intent(in) :: vect(*)
    integer(kind=8), intent(out) :: x1
    integer(kind=8), optional, intent(out) :: x2, x3, x4, x5, x6, x7, x8
    integer(kind=8), optional, intent(out) :: x9, x10, x11, x12, x13, x14
    integer(kind=8), optional, intent(out) :: x15, x16, x17
!-----------------------------------------------------------------------
!
    x1 = vect(1)
    if (present(x2)) x2 = vect(2)
    if (present(x3)) x3 = vect(3)
    if (present(x4)) x4 = vect(4)
    if (present(x5)) x5 = vect(5)
    if (present(x6)) x6 = vect(6)
    if (present(x7)) x7 = vect(7)
    if (present(x8)) x8 = vect(8)
    if (present(x9)) x9 = vect(9)
    if (present(x10)) x10 = vect(10)
    if (present(x11)) x11 = vect(11)
    if (present(x12)) x12 = vect(12)
    if (present(x13)) x13 = vect(13)
    if (present(x14)) x14 = vect(14)
    if (present(x15)) x15 = vect(15)
    if (present(x16)) x16 = vect(16)
    if (present(x17)) x17 = vect(17)
end subroutine
