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

subroutine dplegen(degre, s, l, dlegen)
    implicit none
!
!   DERIVEES DES POLYNOMES DE LEGENDRE
!
#define dpleg2(x)   3.d0*(x)
#define dpleg3(x)   (15.d0*(x)**2-3.d0)/2.d0
#define dpleg4(x)   (35.d0*(x)**3-15.d0*(x))/2.d0
#define dpleg5(x)   15.d0*(21.d0*(x)**4-14.d0*(x)**2+1.d0)/8.d0
#define dpleg6(x)   21.d0*(33.d0*(x)**5-30.d0*(x)**3+5.d0*(x))/8.d0
#define dpleg7(x)   7.d0*(429.d0*(x)**6-495.d0*(x)**4+135.d0*(x)**2-5.d0)/16.d0
!
#include "asterfort/assert.h"
!
    integer(kind=8), intent(in) :: degre
    real(kind=8), intent(in) :: s, l
    real(kind=8), intent(out) :: dlegen
!
!
!     ------------------------------------------------------------------
!
! FONCTION REALISEE:
!
!      CALCUL DES VALEURS DES DERIVEES DES POLYNOMES DE LEGENDRE
!
!     ------------------------------------------------------------------
! ENTREE:
!        DEGRE      : DEGRE DU POLYNOME DE LEGENDRE CONSIDERE
!        S          : ABSCISSE CURVILIGNE DU POINT CONSIDERE
!        L          : LONGUEUR DE LA FISSURE
!
! SORTIE:
!        DLEGEN     : VALEUR DE LA DERIVEE DU POLYNOME DE LEGENDRE
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
        coef = (2.d0/l)*sqrt(1.d0/l)
        dlegen = 0.d0
!
    case (1)
        coef = (2.d0/l)*sqrt(3.d0/l)
        dlegen = coef
!
    case (2)
        coef = (2.d0/l)*sqrt(5.d0/l)
        dlegen = coef*dpleg2(ksi)
!
    case (3)
        coef = (2.d0/l)*sqrt(7.d0/l)
        dlegen = coef*dpleg3(ksi)
!
    case (4)
        coef = (2.d0/l)*sqrt(9.d0/l)
        dlegen = coef*dpleg4(ksi)
!
    case (5)
        coef = (2.d0/l)*sqrt(11.d0/l)
        dlegen = coef*dpleg5(ksi)
!
    case (6)
        coef = (2.d0/l)*sqrt(13.d0/l)
        dlegen = coef*dpleg6(ksi)
!
    case (7)
        coef = (2.d0/l)*sqrt(15.d0/l)
        dlegen = coef*dpleg7(ksi)
!
    case default
        ASSERT(.false.)
!
    end select
!
end subroutine
