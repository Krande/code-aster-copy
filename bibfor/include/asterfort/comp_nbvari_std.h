! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

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
    subroutine comp_nbvari_std(rela_comp, defo_comp, type_cpla, &
                               kit_comp, post_iter, regu_visc, post_incr, &
                               nbVari, numeLaw)
        character(len=16), intent(in) :: rela_comp, defo_comp, type_cpla
        character(len=16), intent(in) :: kit_comp(4), post_iter, regu_visc, post_incr
        integer(kind=8), intent(out) :: nbVari, numeLaw
    end subroutine comp_nbvari_std
end interface
