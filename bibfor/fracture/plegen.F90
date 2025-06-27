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

subroutine plegen(degre, s, l, legen)
    implicit none
!
!   POLYNOMES DE LEGENDRE
!
#define pleg2(x)   (3.d0*(x)*(x)-1.d0)/2.d0
#define pleg3(x)   (x)*(5.d0*(x)*(x)-3.d0)/2.d0
#define pleg4(x)   (35.d0*(x)**4-30.d0*(x)*(x)+3.d0)/8.d0
#define pleg5(x)   (x)*(63.d0*(x)**4-70.d0*(x)*(x)+15.d0)/8.d0
#define pleg6(x)   (231.d0*(x)**6-315.d0*(x)**4+105.d0*(x)*(x)-5.d0)/16.d0
#define pleg7(x)   (x)*(429.d0*(x)**6-693.d0*(x)**4+315.d0*(x)*(x)-35.d0)/16.d0
!
#include "asterfort/assert.h"
!
    integer(kind=8), intent(in) :: degre
    real(kind=8), intent(in) :: s, l
    real(kind=8), intent(out) :: legen
!
!
!     ------------------------------------------------------------------
!
! FONCTION REALISEE:
!
!      CALCUL DES VALEURS DES POLYNOMES DE LEGENDRE
!
!     ------------------------------------------------------------------
! ENTREE:
!        DEGRE      : DEGRE DU POLYNOME DE LEGENDRE CONSIDERE
!        S          : ABSCISSE CURVILIGNE DU POINT CONSIDERE
!        L          : LONGUEUR DE LA FISSURE
!
! SORTIE:
!        LEGEN      : VALEUR DU POLYNOME DE LEGENDRE
!     ------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------

    real(kind=8) :: coef, ksi
!
    ksi = 2.d0*s/l-1.d0
!
    select case (degre)
!
    case (0)
        coef = sqrt(1.d0/l)
        legen = coef*1.d0
!
    case (1)
        coef = sqrt(3.d0/l)
        legen = coef*ksi
!
    case (2)
        coef = sqrt(5.d0/l)
        legen = coef*pleg2(ksi)
!
    case (3)
        coef = sqrt(7.d0/l)
        legen = coef*pleg3(ksi)
!
    case (4)
        coef = sqrt(9.d0/l)
        legen = coef*pleg4(ksi)
!
    case (5)
        coef = sqrt(11.d0/l)
        legen = coef*pleg5(ksi)
!
    case (6)
        coef = sqrt(13.d0/l)
        legen = coef*pleg6(ksi)
!
    case (7)
        coef = sqrt(15.d0/l)
        legen = coef*pleg7(ksi)
!
    case default
        ASSERT(.false.)
!
    end select
!
end subroutine
