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
#include "asterf_types.h"
!
interface
    subroutine greihm(ndim, mecani, press1,&
                      press2, tempe, dimdef, dimcon)
        integer(kind=8) :: ndim
        integer(kind=8) :: mecani(8)
        integer(kind=8) :: press1(9)
        integer(kind=8) :: press2(9)
        integer(kind=8) :: tempe(5)
        integer(kind=8) :: dimdef
        integer(kind=8) :: dimcon
    end subroutine greihm
end interface
