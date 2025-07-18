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
function appmasl(ds_thm, appmasm, &
                 phi, phim, &
                 satur, saturm, &
                 rho, rhom)
!
    use THM_type
!
    implicit none
!
!
    type(THM_DS), intent(in) :: ds_thm
    real(kind=8), intent(in) :: appmasm
    real(kind=8), intent(in) :: phi, phim
    real(kind=8), intent(in) :: satur, saturm
    real(kind=8), intent(in) :: rho, rhom
    real(kind=8) :: appmasl
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Compute quantity of mass from change of volume, porosity and saturation
!
! --------------------------------------------------------------------------------------------------
!
! In  appmasm          : initial quantity of mass
! In  phi              : porosity at end of current time step
! In  phim             : porosity at beginning of current time step
! In  satur            : saturation at end of current time step
! In  saturm           : saturation at beginning of current time step
! In  rho              : volumic mass at end of current time step
! In  rhom             : volumic mass at beginning of current time step
! Out appmasl           : quantity of mass at end of current time step
!
! --------------------------------------------------------------------------------------------------
!
    appmasl = appmasm+phi*satur*rho-phim*saturm*rhom
!
end function
