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
    subroutine dimthm(ds_thm, l_vf, ndim, ndlno, ndlnm)
        use THM_type
        type(THM_DS), intent(in) :: ds_thm
        aster_logical, intent(in) :: l_vf
        integer(kind=8), intent(in)  :: ndim
        integer(kind=8), intent(out) :: ndlno, ndlnm
    end subroutine dimthm
end interface
