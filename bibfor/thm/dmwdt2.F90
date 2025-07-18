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
function dmwdt2(rho11, alp11, phids, satur, cs, &
                dpvpt)
!
    implicit none
!
#include "asterf_types.h"
!
    real(kind=8), intent(in) :: rho11, alp11, phids, satur, cs, dpvpt
    real(kind=8) :: dmwdt2
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Derivative of quantity of mass by liquid pressure - With steam only
!
! --------------------------------------------------------------------------------------------------
!
! In  rho11            : volumic mass of liquid
! In  alp11            : thermic dilatation of liquid
! In  phids            : product porosity by derivative of saturation (/pc)
! In  satur            : saturation
! In  cs               : Biot modulus of solid matrix
! In  dpvpt            : derivative of steam pressure by temperature
! Out dmwdt2           : derivative of quantity of mass by liquid pressure
!
! --------------------------------------------------------------------------------------------------
!
    dmwdt2 = rho11*(-3.d0*alp11+(phids-satur*satur*cs+satur*cs)*dpvpt)
!
end function
