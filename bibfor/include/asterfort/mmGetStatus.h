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
    subroutine mmGetStatus(indco       ,&
                           l_prev_cont_, l_prev_fric_ ,&
                           indco_prev_ , indadhe_prev_, indadhe2_prev_)
        integer(kind=8), intent(out) :: indco
        aster_logical, optional, intent(out) :: l_prev_cont_, l_prev_fric_
        integer(kind=8), optional, intent(out) :: indco_prev_, indadhe_prev_, indadhe2_prev_
    end subroutine mmGetStatus
end interface
