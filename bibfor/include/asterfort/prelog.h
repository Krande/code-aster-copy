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
    subroutine prelog(ndim, lgpg, vim, gn, lamb,&
                      logl, fPrev, fCurr, epslPrev, epslIncr,&
                      tlogPrev, lCorr, iret)
        integer(kind=8), intent(in) :: ndim, lgpg
        real(kind=8), intent(in) :: vim(lgpg)
        real(kind=8), intent(in) :: fPrev(3, 3), fCurr(3, 3)
        real(kind=8), intent(out) :: epslPrev(6), epslIncr(6)
        real(kind=8), intent(out) :: tlogPrev(6)
        real(kind=8), intent(out) :: gn(3, 3), lamb(3), logl(3)
        aster_logical, intent(in) :: lCorr
        integer(kind=8), intent(out) :: iret
    end subroutine prelog
end interface
