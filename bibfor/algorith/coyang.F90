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

function coyang(dist, dteta, rayon, omega, uc, &
                uct, l, lt)
! CALCUL DE LA FONCTION DE COHERENCE SUIVANT AU-YANG
    implicit none
!
! *****************   DECLARATIONS DES VARIABLES   ********************
!
!
! ARGUMENTS
! ---------
#include "jeveux.h"
#include "asterc/r8pi.h"
    real(kind=8) :: dist, dteta, rayon, omega, uc, uct, l, lt, dteta2, pi
!
!
! VARIABLES LOCALES
! -----------------
    real(kind=8) :: coyang, codist, coteta
!
! si la vitesse convective transversale est supérieure a 10 fois la vitesse convective
! longitudinale alors la convection des tourbillons transversale n est pas prise
! en compte
! CALCUL DE LA FONCTION DE COHERENCE
    pi = r8pi()
    dteta2 = min(abs(dteta), abs(dteta+2*pi), abs(dteta-2*pi))
    if (uct .gt. 10*uc) then
        codist = exp(-abs(dist)/l)*cos(omega*dist/uc)
        coteta = exp(-rayon*dteta2/lt)
        coyang = codist*coteta
    else
        codist = exp(-abs(dist)/l)*exp(-rayon*abs(dteta2)/lt)
        coteta = cos(rayon*omega*dteta/uct+omega*dist/uc)
        coyang = codist*coteta
    end if
!
end function
