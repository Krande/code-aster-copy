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
subroutine sigmap(ds_thm, &
                  satur, signe, tbiot, dp2, dp1, dpi, &
                  sigmp)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: signe, tbiot(6), satur, dp1, dp2, dpi
    real(kind=8), intent(out) :: sigmp(6)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute mechanical stresses from pressures
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  signe            : sign for saturation
! In  tbiot            : tensor of Biot
! In  satur            : value of saturation
! In  dp1              : increment of capillary pressure
! In  dp2              : increment of gaz pressure
! In  dpi              : variation of the hydraulic pressure at end of current time
! Out sigmp            : stress component for pressure
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i
!
! --------------------------------------------------------------------------------------------------
!
    do i = 1, 6
        if (ds_thm%ds_behaviour%l_stress_bishop) then
            if ((ds_thm%ds_behaviour%rela_hydr) .eq. 'HYDR_TABBAL') then
                sigmp(i) = -tbiot(i)*dpi

            else
                sigmp(i) = tbiot(i)*satur*signe*dp1-tbiot(i)*dp2
            end if
        else
            sigmp(i) = -tbiot(i)*dp2
        end if
    end do
!
end subroutine
