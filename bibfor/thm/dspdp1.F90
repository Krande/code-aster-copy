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
subroutine dspdp1(ds_thm, signe, tbiot, satur, dsdp1, phi0, ep, surf, sbjh, wbjh)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: signe, tbiot(6), satur
    real(kind=8), intent(out) :: dsdp1(6)
    real(kind=8), intent(in)  :: phi0, ep, surf, sbjh, wbjh
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Derivative of _pressure part_ of stresses by capillary pressure
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  signe            : sign for saturation
! In  tbiot            : tensor of Biot
! In  satur            : value of saturation
! Out dsdp1            : derivative of pressure part of stress by capillary pressure
! In  phi0             : initial porosity (THM_INIT)
! In  ep               : thickness of the adsorbed water layer
! In  surf             : specific surface of the material
! In  sbjh             : saturated pores volume fraction from BJH - At end of current step
! In  wbjh             : unsaturatred pores surface fraction from BJH - At end of current step
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i
!
! --------------------------------------------------------------------------------------------------
!
    do i = 1, 6
        if (ds_thm%ds_behaviour%l_stress_bishop) then
            if ((ds_thm%ds_behaviour%rela_hydr) .eq. 'HYDR_TABBAL') then
                dsdp1(i) = -tbiot(i)*(-sbjh-((2./3.)*(surf/phi0)*wbjh*ep))

            else
                dsdp1(i) = signe*tbiot(i)*satur
            end if
        else
            dsdp1(i) = 0.d0
        end if
    end do
!
end subroutine
