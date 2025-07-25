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
subroutine zerop3(a, b, c, x, n)
!
!
! aslint: disable=
    implicit none
#include "asterc/r8pi.h"
    real(kind=8) :: a, b, c, x(3)
    integer(kind=8) :: n
!
! ----------------------------------------------------------------------
! RESOLUTION D'UN POLYNOME DE DEGRE 3 : X**3 + A X**2 + B X + C = 0
! ----------------------------------------------------------------------
! IN  A,B,C COEFFICIENTS DU POLYNOME
! OUT X     RACINES DANS L'ORDRE DECROISSANT
! OUT N     NOMBRE DE RACINES
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i
    real(kind=8) :: p, q, delta
    real(kind=8) :: tau, cs, alpha, u, t, y(3)
    real(kind=8) :: pi, v1, v2
!
! -- ON SE RAMENE A : Y**3 - P Y - Q = 0   AVEC Y = X + A/3
!
    p = a**2/3.d0-b
    q = a*b/3-c-2*a**3/27.d0
!
! -- TRAITEMENT DES CAS PARTICULIERS
!
    if (p .eq. 0 .and. q .eq. 0) then
        n = 3
        y(1) = 0
        y(2) = 0
        y(3) = 0
        goto 1000
    else if (p .eq. 0) then
        n = 1
        v1 = abs(q)**(1.d0/3.d0)
        y(1) = sign(v1, q)
        goto 1000
    else if (p .lt. 0 .and. q .eq. 0) then
        n = 1
        y(1) = 0
        goto 1000
    else if (p .gt. 0 .and. q .eq. 0) then
        n = 3
        y(1) = sqrt(p)
        y(2) = 0
        y(3) = -sqrt(p)
        goto 1000
    end if
!
! -- SOLUTION UNIQUE SI P<0  OU  ABS(Q) > 2 (P/3) ** 3/2
!
    if (p .lt. 0 .or. abs(q) .gt. 2*abs(p/3)**1.5d0) then
        n = 1
        delta = 27*q**2-4*p**3
        t = (27*q+sign(sqrt(abs(27*delta)), q))/2
        v2 = abs(t)**(1.d0/3.d0)
        u = sign(v2, t)
        y(1) = p/u+u/3
!
! -- SINON : TROIS RACINES
!
    else
        n = 3
        pi = r8pi()
        tau = 2*sqrt(p/3.d0)
        cs = 4*q/tau**3
        if (cs .ge. 1) then
            alpha = 0
        else if (cs .le. -1) then
            alpha = pi/3
        else
            alpha = atan2(sqrt(1.d0-(cs**2.d0)), cs)
            alpha = alpha/3.d0
        end if
!
        if (alpha .le. pi/3) then
            y(1) = tau*cos(alpha)
            y(2) = tau*cos(alpha-2*pi/3.d0)
            y(3) = tau*cos(alpha+2*pi/3.d0)
        else if (alpha .le. 2*pi/3.d0) then
            y(1) = tau*cos(alpha-2*pi/3.d0)
            y(2) = tau*cos(alpha)
            y(3) = tau*cos(alpha+2*pi/3.d0)
        else
            y(1) = tau*cos(alpha-2*pi/3.d0)
            y(2) = tau*cos(alpha+2*pi/3.d0)
            y(3) = tau*cos(alpha)
        end if
    end if
!
1000 continue
!
    do i = 1, n
        x(i) = y(i)-a/3
    end do
!
end subroutine
