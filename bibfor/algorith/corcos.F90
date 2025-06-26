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

function corcos(d1, d2, mes1, mes2, uc, uct, l, &
                lt, omega)
    implicit none
!corcos(d1,d2,local(1), local(2), uc, uct, l,lt,omega)
! *****************   DECLARATIONS DES VARIABLES   ********************
!
!
! ARGUMENTS
! ---------
#include "jeveux.h"
    real(kind=8) :: d1, d2, omega, l, lt, uc, uct, mes1, mes2
!
!
! VARIABLES LOCALES
! -----------------
    real(kind=8) :: corcos
!
!
!
! CALCUL DE LA FONCTION DE COHERENCE
! Attention : Si la vitesse convective transversale est supérieure a 10 fois
! la vitesse convective longitudinale alors la convection des tourbillons
! transversale n est pas prise en compte
!
    if (uct .gt. 10*uc) then
        corcos = exp(-d2/lt)*exp(-d1/l)*cos(omega*mes1/uc)
    else
        corcos = exp(-d2/lt)*exp(-d1/l)*cos(omega*mes1/uc+omega*mes2/uct)
    end if
!
end function
