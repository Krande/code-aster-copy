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
!
#include "asterf_types.h"
!
interface
    subroutine iniMat0(l, x1, x2, x3, x4, x5,&
                   x6, x7, x8, x9, x10, x11, x12, x13, x14,&
                   x15, x16, x17, x18, x19)
        integer(kind=8), intent(in) :: l
        real(kind=8), intent(out) :: x1(l,l)
        real(kind=8), optional, intent(out) :: x2(l,l)
        real(kind=8), optional, intent(out) :: x3(l,l)
        real(kind=8), optional, intent(out) :: x4(l,l)
        real(kind=8), optional, intent(out) :: x5(l,l)
        real(kind=8), optional, intent(out) :: x6(l,l)
        real(kind=8), optional, intent(out) :: x7(l,l)
        real(kind=8), optional, intent(out) :: x8(l,l)
        real(kind=8), optional, intent(out) :: x9(l,l)
        real(kind=8), optional, intent(out) :: x10(l,l)
        real(kind=8), optional, intent(out) :: x11(l,l)
        real(kind=8), optional, intent(out) :: x12(l,l)
        real(kind=8), optional, intent(out) :: x13(l,l)
        real(kind=8), optional, intent(out) :: x14(l,l)
        real(kind=8), optional, intent(out) :: x15(l,l)
        real(kind=8), optional, intent(out) :: x16(l,l)
        real(kind=8), optional, intent(out) :: x17(l,l)
        real(kind=8), optional, intent(out) :: x18(l,l)
        real(kind=8), optional, intent(out) :: x19(l,l)
    end subroutine iniMat0
end interface
