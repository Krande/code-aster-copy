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
subroutine pj3da3(m, a, b, c, ok, &
                  la, lb, lc, d2)
    implicit none
#include "asterf_types.h"
    real(kind=8) :: m(3), a(3), b(3), c(3), d2, la, lb, lc
    aster_logical :: ok
! BUT :
!   TROUVER LES COORDONNEES BARYCENTRIQUES (LA,LB,LC) DE LA PROJECTION P
!   D'UN POINT M SUR UN TRIANGLE (A,B,C) .
!
!  IN   M(3)    R : COORDONNEES DE M
!  IN   A(3)    R : COORDONNEES DE A
!  IN   B(3)    R : COORDONNEES DE B
!  IN   C(3)    R : COORDONNEES DE C
!
!  OUT  OK         L  :/.TRUE.   : P EST INTERIEUR AU TRIANGLE ABC
!                      /.FALSE.  : P EST EXTERIEUR AU TRIANGLE
!                       (SI .FALSE.   D2 N'EST PAS CALCULE)
!  OUT  D2         R  : CARRE DE LA DISTANCE ENTRE M ET P
!  OUT  LA,LB,LC   R  : COORDONNEES BARYCENTRIQUES DE P DANS ABC
!
!
! ----------------------------------------------------------------------
    integer(kind=8) :: k
    real(kind=8) :: delta, p(3)
    real(kind=8) :: ab(3), ac(3), am(3), a11, a22, a12, b1, b2
! DEB ------------------------------------------------------------------
    do k = 1, 3
        ab(k) = b(k)-a(k)
        ac(k) = c(k)-a(k)
        am(k) = m(k)-a(k)
    end do
!
    a11 = ab(1)*ab(1)+ab(2)*ab(2)+ab(3)*ab(3)
    a22 = ac(1)*ac(1)+ac(2)*ac(2)+ac(3)*ac(3)
    a12 = ab(1)*ac(1)+ab(2)*ac(2)+ab(3)*ac(3)
!
    b1 = ab(1)*am(1)+ab(2)*am(2)+ab(3)*am(3)
    b2 = ac(1)*am(1)+ac(2)*am(2)+ac(3)*am(3)
!
    delta = a11*a22-a12*a12
    if (delta .eq. 0.d0) then
        ok = .false.
        goto 999
    end if
    lb = (a22*b1-a12*b2)/delta
    lc = (a11*b2-a12*b1)/delta
    la = 1.d0-lb-lc
!
    if ((la .ge. 0.d0) .and. (la .le. 1.d0) .and. (lb .ge. 0.d0) .and. (lb .le. 1.d0) .and. &
        (lc .ge. 0.d0) .and. (lc .le. 1.d0)) then
        ok = .true.
        do k = 1, 3
            p(k) = la*a(k)+lb*b(k)+lc*c(k)
            p(k) = m(k)-p(k)
        end do
        d2 = p(1)*p(1)+p(2)*p(2)+p(3)*p(3)
    else
        ok = .false.
    end if
!
999 continue
end subroutine
