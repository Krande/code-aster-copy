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
    subroutine nd_mstp_time(ds_inout, list_func_acti, time_prev_step, l_comp_mstp)
        use NonLin_Datastructure_type
        type(NL_DS_InOut), intent(in) :: ds_inout
        integer(kind=8), intent(in) :: list_func_acti(*)
        real(kind=8), intent(out) :: time_prev_step
        aster_logical, intent(out) :: l_comp_mstp
    end subroutine nd_mstp_time
end interface
