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
#include "asterf_types.h"
!
interface
    subroutine xcoocy(ndim, xg, pfon, p, rg, tg, l_not_zero)
        integer(kind=8) :: ndim
        real(kind=8) :: rg
        real(kind=8) :: tg
        real(kind=8) :: xg(ndim)
        real(kind=8) :: pfon(ndim)
        real(kind=8) :: p(3,3)
        aster_logical :: l_not_zero
    end subroutine xcoocy
end interface
