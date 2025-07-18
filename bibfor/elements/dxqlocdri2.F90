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

subroutine dxqlocdri2(btgmemb, matloc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
    real(kind=8) :: btgmemb(*)
    real(kind=8) :: matloc(*)
!----------------------------------------------------------
!     IN  GMEMB  : MATRICE DE MEMBRANE (DRILLING) CARREE
!     OUT MATLOC : MATRICE DE RIGIDITE OU DE MASSE LOCALE
!                  REMPLISSAGE DE MATELEM LOCAL (300 TERMES) AVEC
!                  LES TERMES DE MEMBRANE DRILLING DRZ
!----------------------------------------------------------
!
!
    integer(kind=8) :: im(32), jm(32)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: k
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------

    data jm/&
     &   16, 17, 27, 34, 67, 68, 73, 74, 84, 90,&
     &   97, 103, 154, 155, 160, 161, 166, 167, 177, 183,&
     &  189, 196, 202, 208, 277, 278, 283, 284, 289, 290,&
     &  295, 296/

    data im/&
     &    1, 2, 3, 4, 9, 10, 11, 12, 5, 13,&
     &    6, 14, 17, 18, 19, 20, 21, 22, 7, 15,&
     &   23, 8, 16, 24, 25, 26, 27, 28, 29, 30,&
     &   31, 32/

!     ------------------------------------------------------------------

!                       ---- TERMES DE MEMBRANE (DRILL OUT-OF-DIAGONAL)
    do k = 1, 32
        matloc(jm(k)) = matloc(jm(k))+btgmemb(im(k))
    end do

end subroutine
