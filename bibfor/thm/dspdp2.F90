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
subroutine dspdp2(ds_thm, tbiot, dsdp2)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: tbiot(6)
    real(kind=8), intent(out) :: dsdp2(6)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Derivative of pressure part of stress by gaz pressure
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  tbiot            : tensor of Biot
! Out dsdp2            : derivative of pressure part of stress by gaz pressure
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i
!
! --------------------------------------------------------------------------------------------------
!
    do i = 1, 6
        if (ds_thm%ds_behaviour%l_stress_bishop) then
            dsdp2(i) = -tbiot(i)
        else
            dsdp2(i) = -tbiot(i)
        end if
    end do
!
end subroutine
