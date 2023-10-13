! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine addVecLumped(vec, vec_sub, ise, size, connec)
!
    implicit none
!
#include "FE_module.h"
!
    real(kind=8), intent(inout) :: vec(MAX_BS)
    real(kind=8), intent(in) :: vec_sub(MAX_BS)
    integer, intent(in) :: ise, size, connec(4,27)
!
! --------------------------------------------------------------------------------------------------
!
! Add subCell lumped vector to global vector
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, idi
!
! --------------------------------------------------------------------------------------------------
!
    do i = 1, size
        idi = connec(ise, i)
        vec(idi) = vec(idi)+vec_sub(i)
    end do
!
end subroutine
