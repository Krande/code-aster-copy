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

subroutine dxqlocdri1(gmemb, matloc)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
    real(kind=8) :: gmemb(*)
    real(kind=8) :: matloc(*)
!----------------------------------------------------------
!     IN  GMEMB  : MATRICE DE MEMBRANE (DRILLING) CARREE
!     OUT MATLOC : MATRICE DE RIGIDITE OU DE MASSE LOCALE
!                  REMPLISSAGE DE MATELEM LOCAL (300 TERMES) AVEC
!                  LES TERMES DE MEMBRANE DRILLING DRZ
!----------------------------------------------------------
!
    integer(kind=8) :: im(10), jm(10)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: k
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
    data jm/&
     &   21, 72, 78, 159, 165, 171, 282, 288, 294, 300/
    data im/&
     &    1, 5, 6, 9, 10, 11, 13, 14, 15, 16/

!     ------------------------------------------------------------------

!                          ---- TERMES DE MEMBRANE (DRILLING-DIAGONAL)

    do k = 1, 10
        matloc(jm(k)) = matloc(jm(k))+gmemb(im(k))
    end do

end subroutine
