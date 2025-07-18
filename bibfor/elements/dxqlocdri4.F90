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

subroutine dxqlocdri4(bxb, matloc)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
    real(kind=8) :: bxb(*)
    real(kind=8) :: matloc(*)
!----------------------------------------------------------
!     IN  GMEMB  : MATRICE DE MEMBRANE (DRILLING) CARREE
!     OUT MATLOC : MATRICE DE RIGIDITE OU DE MASSE LOCALE
!                  REMPLISSAGE DE MATELEM LOCAL (300 TERMES) AVEC
!                  LES TERMES DE MEMBRANE DRILLING DRZ
!----------------------------------------------------------
!
    integer(kind=8) :: im(78), jm(78)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: k
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------

    data jm/&
     &    1, 2, 3, 16, 17, 21, 22, 23, 27, 28,&
     &   29, 30, 34, 35, 36, 67, 68, 72, 73, 74,&
     &   78, 79, 80, 84, 85, 86, 90, 91, 92, 93,&
     &   97, 98, 99, 103, 104, 105, 154, 155, 159, 160,&
     &  161, 165, 166, 167, 171, 172, 173, 177, 178, 179,&
     &  183, 184, 185, 189, 190, 191, 192, 196, 197, 198,&
     &  202, 203, 204, 208, 209, 210, 277, 278, 282, 283,&
     &  284, 288, 289, 290, 294, 295, 296, 300/

    data im/&
     &    1, 13, 14, 25, 26, 27, 37, 38, 39, 40,&
     &   49, 50, 51, 52, 53, 61, 62, 63, 64, 65,&
     &   66, 73, 74, 75, 76, 77, 78, 79, 85, 86,&
     &   87, 88, 89, 90, 91, 92, 97, 98, 99, 100,&
     &  101, 102, 103, 104, 105, 109, 110, 111, 112, 113,&
     &  114, 115, 116, 117, 118, 121, 122, 123, 124, 125,&
     &  126, 127, 128, 129, 130, 131, 133, 134, 135, 136,&
     &  137, 138, 139, 140, 141, 142, 143, 144/

!     ------------------------------------------------------------------

!                       ---- TERMES DE MEMBRANE
    do k = 1, 78
        matloc(jm(k)) = matloc(jm(k))+bxb(im(k))
    end do

end subroutine
