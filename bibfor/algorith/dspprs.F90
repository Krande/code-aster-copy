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

function dspprs(k, u, d, rho, f, &
                fcoup)
    implicit none
!
! *****************   DECLARATIONS DES VARIABLES   ********************
!
!
! ARGUMENTS
! ---------
#include "jeveux.h"
    real(kind=8) :: k, u, d, rho, f, dspprs, fcoup
!
!
! VARIABLES LOCALES
! -----------------
    real(kind=8) :: dsp, fr
!
! CALCUL DE LA FREQUENCE REDUITE
    fr = f*d/u
    if (fr .le. fcoup) then
        dsp = (k**2)*((rho*(u**2))**2)*(d**3)
    else
        dsp = 0.d0
    end if
    dspprs = dsp
!
end function
