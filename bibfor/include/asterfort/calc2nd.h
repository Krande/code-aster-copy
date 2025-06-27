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
    subroutine calc2nd( ds_thm, j_mater, &
                        lMatr, lSigm, &
                        ndim, dimdef, dimcon, &
                        adde2nd, adco2nd, &
                        defgem, defgep, &
                        congem, congep, &
                        dsde)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        aster_logical, intent(in) :: lMatr, lSigm
        integer(kind=8), intent(in) :: j_mater, ndim, dimdef, dimcon
        integer(kind=8), intent(in) :: adde2nd, adco2nd
        real(kind=8), intent(in) :: defgem(dimdef)
        real(kind=8), intent(in) :: defgep(dimdef)
        real(kind=8), intent(in) :: congem(dimcon)
        real(kind=8), intent(inout) :: congep(dimcon)
        real(kind=8), intent(inout) :: dsde(dimcon, dimdef)
    end subroutine calc2nd
end interface 
