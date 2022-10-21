! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
    subroutine setMat33Tab(tab, x1, x2, x3, x4, x5,&
                      x6, x7, x8, x9, x10, x11, x12, x13, x14,&
                      x15, x16, x17)
        real(kind=8), intent(inout) :: tab(3,3,*)
        real(kind=8), intent(in) :: x1(3,3)
        real(kind=8), optional, intent(in) :: x2(3,3)
        real(kind=8), optional, intent(in) :: x3(3,3)
        real(kind=8), optional, intent(in) :: x4(3,3)
        real(kind=8), optional, intent(in) :: x5(3,3)
        real(kind=8), optional, intent(in) :: x6(3,3)
        real(kind=8), optional, intent(in) :: x7(3,3)
        real(kind=8), optional, intent(in) :: x8(3,3)
        real(kind=8), optional, intent(in) :: x9(3,3)
        real(kind=8), optional, intent(in) :: x10(3,3)
        real(kind=8), optional, intent(in) :: x11(3,3)
        real(kind=8), optional, intent(in) :: x12(3,3)
        real(kind=8), optional, intent(in) :: x13(3,3)
        real(kind=8), optional, intent(in) :: x14(3,3)
        real(kind=8), optional, intent(in) :: x15(3,3)
        real(kind=8), optional, intent(in) :: x16(3,3)
        real(kind=8), optional, intent(in) :: x17(3,3)
    end subroutine setMat33Tab
end interface
