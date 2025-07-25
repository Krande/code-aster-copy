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
subroutine zerofb(func, x1, x2, tol, itmax, &
                  zbrent, iret, iter)
!
    implicit none
!
#include "asterc/r8prem.h"
    interface
        function func(x)
            real(kind=8) :: func, x
        end function func
    end interface
    integer(kind=8) :: itmax, iter, iret
    real(kind=8) :: zbrent, tol, x1, x2, eps
    real(kind=8) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
! ----------------------------------------------------------------------
!     BUT : TROUVER LE ZERO D'UNE FONCTION SCALAIRE REELLE
!     AVEC LA METHODE DE BRENT
!
!     USING BRENT'S METHOD, FIND THE ROOT OF A FUNCTION func KNOWN TO
!     LIE BETWEEN X1 AND X2. THE ROOT, RETURNED AS ZBRENT, WILL BE
!     REFINED UNTIL ITS ACCURACY IS TOL.
!     PARAMETERS: MAXIMUM ALLOWED NUMBER OF ITERATIONS
!
! IN  func       : FONCTION func
! IN  X1, X2  : INTERVELLE DE RECHERCHE
! IN  TOL     : PRECISION ABSOLUE : LA SOLUTION X EST TELLE QUE func(X)<TOL
! IN  ITMAX   : NOMBRE D'ITERATIONS MAXIMUM
! OUT ZBRENT  : ZERO DE func
! OUT IRET    : CODE RETOUR : IRET = 0 : OK
!             :               IRET = 1 : NITER INSUFFISANT OU AUTRE PB
! OUT ITER    : NOMBRE D'ITERATIONS EFFECTUEES
! ----------------------------------------------------------------------
!
    eps = r8prem()
    iret = 0
    iter = 0
    a = x1
    b = x2
    fa = func(a)
    fb = func(b)
!
    if (fa .gt. 0.d0 .and. fb .gt. 0.d0 .or. fa .lt. 0.d0 .and. fb .lt. 0.d0) then
!
        iret = 1
        goto 999
!
    end if
!
    c = b
    fc = fb
!
    do iter = 1, itmax
!
        if (fb .gt. 0.d0 .and. fc .gt. 0.d0 .or. fb .lt. 0.d0 .and. fc .lt. 0.d0) then
!         RENAME A, B, C AND ADJUST BOUNDING INTERVAL D.
            c = a
            fc = fa
            d = b-a
            e = d
        end if
!
        if (abs(fc) .lt. abs(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end if
!
!       CONVERGENCE CHECK.
        tol1 = 2.d0*eps*abs(b)
        xm = 0.5d0*(c-b)
        if (abs(xm) .le. tol1 .or. abs(fb) .lt. tol) then
            zbrent = b
            goto 999
        end if
!
        if (abs(e) .ge. tol1 .and. abs(fa) .gt. abs(fb)) then
!         ATTEMPT INVERSE QUADRATIC INTERPOLATION.
            s = fb/fa
            if (a .eq. c) then
                p = 2.d0*xm*s
                q = 1.d0-s
            else
                q = fa/fc
                r = fb/fc
                p = s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
                q = (q-1.d0)*(r-1.d0)*(s-1.d0)
            end if
!         CHECK WHETHER IN BOUNDS.
            if (p .gt. 0.d0) q = -q
            p = abs(p)
            if (2.d0*p .lt. min(3.d0*xm*q-abs(tol1*q), abs(e*q))) then
!           ACCEPT INTERPOLATION.
                e = d
                d = p/q
            else
!           INTERPOLATION FAILED, USE BISECTION.
                d = xm
                e = d
            end if
        else
!         BOUNDS DECREASING TOO SLOWLY, USE BISECTION.
            d = xm
            e = d
        end if
!
!       MOVE LAST BEST GUESS TO A.
        a = b
        fa = fb
!
!       EVALUATE NEW TRIAL ROOT.
        if (abs(d) .gt. tol1) then
            b = b+d
        else
            b = b+sign(tol1, xm)
        end if
        fb = func(b)
    end do
!
    iret = 1
    zbrent = b
!
999 continue
end subroutine
