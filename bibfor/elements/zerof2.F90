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
subroutine zerof2(func, x0, xap, epsi, nitmax, &
                  solu, iret, n)
    implicit none
#include "asterfort/utmess.h"
!
!     ARGUMENTS:
!     ----------
    interface
        function func(x)
            real(kind=8) :: func, x
        end function func
    end interface
    real(kind=8) :: x0, xap, epsi, solu
    integer(kind=8) :: nitmax, iret
! ----------------------------------------------------------------------
!     BUT:
!         TROUVER UNE RACINE DE L'EQUATION func(X)=0
!         ON SUPPOSE QUE LA FONCTION func EST CROISSANTE ET QUE func(X0)<0
!         ON EMPLOIE LA METHODE DE SECANTE UTILISEE DANS ZEROFO AVEC
!          EN PLUS UN "COUP" DE DICHOTOMIE TOUS LES 3 ITERATIONS
!          POUR FACILITER LA CONVERGENCE SI func EST TRES NON-LINEAIRE
!
!     IN:
!         func  : FONCTION DONT ON CHERCHE LE "ZERO"
!         X0 : POINT 0
!         XAP: APPROXIMATION DE LA SOLUTION.
!        EPSI: TOLERANCE ABSOLU SUR LE ZERO CHERCHE : ABS(func(SOLU))<EPSI
!      NITMAX: NOMBRE MAXI D'ITERATIONS AUTORISEES.
!
!     OUT:
!         SOLU: VALEUR DE LA RACINE CHERCHEE.
!         IRET: CODE RETOUR DE LA RECHERCHE DE ZERO DE func(X)=0
!                   IRET=0 => PAS DE PROBLEME
!                   IRET=1 => ECHEC
!     N       : NOMBRE D'ITERATIONS REALISEES
! ----------------------------------------------------------------------
    real(kind=8) :: fy, fz, x, y, z, a, b, fa, fb, fdbg(20), xdbg(20), ecresd
    real(kind=8) :: fx
    real(kind=8) :: valr(44)
    integer(kind=8) :: n, k, nd
    integer(kind=8) :: vali
! DEB-------------------------------------------------------------------
!
!     INITIALISATIONS
!
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    n = 1
    x = x0
    fx = func(x0)
    y = xap
    fy = func(y)
!
    if (abs(fy) .lt. epsi) then
        z = y
        goto 90
    end if
!
    if (abs(x-y) .le. 1d-15) then
        goto 100
    end if
!
!     DEBUT DES ITERATIONS
!
10  continue
    if (fy .gt. 0.d0) then
        a = x
        b = y
        fa = fx
        fb = fy
!       ND = INT(SQRT(DBLE(NITMAX)))
        nd = 3
20      continue
        if ((n-(n/nd)*nd) .eq. 0) then
            z = (a+b)*0.5d0
        else
            z = (a*fb-b*fa)/(fb-fa)
        end if
!
        n = n+1
        fz = func(z)
!
        if (abs(fz) .lt. epsi) goto 90
        ecresd = abs(b-a)
! SOLUTION PROVISOIRE PERMETTANT DE PASSER LES CAS
! DIFFICILES CF AL98-193 AL98-197
! IL FAUDRAIT FAIRE MIEUX....
        if (ecresd .le. (epsi*b)) goto 90
        if (n .gt. nitmax) goto 98
        if (fz .lt. 0.d0) then
            a = z
            fa = fz
        else
            b = z
            fb = fz
        end if
        goto 20
    else
!
        if (fy .lt. fx) goto 99
!
        if (fy .eq. fx) then
            goto 100
        end if
!
        z = (x*fy-y*fx)/(fy-fx)
!
!
        if (abs(z-y) .le. 1d-15) then
            goto 100
        end if
!
        n = n+1
        x = y
        fx = fy
        y = z
        fy = func(z)
!
!
        if (abs(fy) .lt. epsi) goto 90
        if (n .gt. nitmax) goto 98
    end if
    goto 10
!
90  continue
    solu = z
    goto 999
!
98  continue
    iret = 1
    goto 999
!
99  continue
    do k = 1, 20
        xdbg(k) = xap/(21-k)
        fdbg(k) = func((xap)/(21-k))
    end do
    vali = n
    valr(1) = x
    valr(2) = fx
    valr(3) = y
    valr(4) = fy
    do k = 1, 20
        valr(4+k) = xdbg(k)
        valr(24+k) = fdbg(k)
    end do
!
    call utmess('F', 'ELEMENTS5_39', si=vali, nr=44, valr=valr)
!
100 continue
    do k = 1, 20
        xdbg(k) = xap/(21-k)
        fdbg(k) = func((xap)/(21-k))
    end do
    vali = n
    valr(1) = x
    valr(2) = fx
    valr(3) = y
    valr(4) = fy
    do k = 1, 20
        valr(4+k) = xdbg(k)
        valr(24+k) = fdbg(k)
    end do
    call utmess('F', 'ELEMENTS5_40', si=vali, nr=44, valr=valr)
!
999 continue
end subroutine
