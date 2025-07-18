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

subroutine zerofo(func, x0, xap, epsi, nitmax, &
                  solu, iret, n)
    implicit none
!
    interface
        function func(x)
            real(kind=8) :: x
            real(kind=8) :: func
        end function
    end interface
    real(kind=8) :: x0, xap, epsi, solu
    integer(kind=8) :: nitmax, iret, n
!
!
! ----------------------------------------------------------------------
!     BUT:
!         TROUVER UNE RACINE DE L'EQUATION func(X)=0
!         ON SUPPOSE QUE LA FONCTION func EST CROISSANTE ET QUE func(X0)<0
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
!     IRET    : CODE RETOUR DE LA RECHERCHE DE ZERO
!               IRET=0 => PAS DE PROBLEME
!               IRET=1 => ECHEC DANS LA RECHERCHE DE ZERO
!     N       : NOMBRE D'ITERATIONS REALISEES
!
! ----------------------------------------------------------------------
    real(kind=8) :: fx, fy, fz, x, y, z, a, b
! DEB-------------------------------------------------------------------
!
!     INITIALISATIONS
!
    iret = 1
    n = 1
    x = x0
    fx = func(x0)
    if (abs(fx) .lt. epsi) then
        z = 0.d0
        goto 800
    end if
    y = xap
    fy = func(y)
!
!     DEBUT DES ITERATIONS
!
10  continue
    if (fy .gt. 0.d0) then
        a = x
        b = y
20      continue
        if (fx .eq. fy) goto 999
        z = y-(y-x)*fy/(fy-fx)
        if (((z-a)*(z-b)) .gt. 0.d0) then
            z = (a+b)/2.d0
        end if
!
        n = n+1
        fz = func(z)
        if (abs(fz) .lt. epsi) goto 800
        if (n .gt. nitmax) goto 999
        if (fz .lt. 0.d0) then
            a = z
        else
            b = z
        end if
        x = y
        fx = fy
        y = z
        fy = fz
        goto 20
    else
        if (fy .lt. fx) goto 999
        if (fx .eq. fy) goto 999
        z = y-(y-x)*fy/(fy-fx)
        n = n+1
        x = y
        fx = fy
        y = z
        fy = func(z)
!
        if (abs(fy) .lt. epsi) goto 800
        if (n .gt. nitmax) goto 999
    end if
    goto 10
!
!     SUCCES
800 continue
    solu = z
    iret = 0
!
!     SORTIE
999 continue
end subroutine
