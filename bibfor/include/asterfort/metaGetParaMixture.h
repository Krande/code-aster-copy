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
    subroutine metaGetParaMixture(poum, fami, kpg, ksp, jvMaterCode, &
                                  l_visc, metaType, nbPhase, zalpha, fmix, &
                                  sy)
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
        integer(kind=8), intent(in) :: metaType
        integer(kind=8), intent(in) :: nbPhase
        aster_logical, intent(in) :: l_visc
        real(kind=8), intent(in) :: zalpha
        real(kind=8), intent(out) :: fmix
        real(kind=8), optional, intent(out) :: sy(nbPhase)
    end subroutine metaGetParaMixture
end interface
