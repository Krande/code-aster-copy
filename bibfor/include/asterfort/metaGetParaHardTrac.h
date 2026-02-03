! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
    subroutine metaGetParaHardTrac(jvMaterCode, metaType, nbPhase, &
                                   l_temp, temp, &
                                   epseq, h0, rp_, nbValeMaxi_)
        integer(kind=8), intent(in) :: jvMaterCode
        integer(kind=8), intent(in) :: metaType
        integer(kind=8), intent(in) :: nbPhase
        aster_logical, intent(in) :: l_temp
        real(kind=8), intent(in) :: temp
        real(kind=8), intent(in) :: epseq(nbPhase)
        real(kind=8), intent(out) :: h0(nbPhase)
        real(kind=8), optional, intent(out) :: rp_(nbPhase)
        integer(kind=8), optional, intent(out) :: nbValeMaxi_
    end subroutine metaGetParaHardTrac
end interface
