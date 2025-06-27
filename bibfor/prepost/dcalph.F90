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
subroutine dcalph(x, y, nbpts, pe)
!
!**********************************************************
!              BUT DE CETTE ROUTINE :                     *
! CALCULER L ORDRE DE LA SINGULARITE PAR LA METHODE       *
! DES MOINDRES CARRES                                     *
! RECHERCHE DE LA PENTE ET DE LA CONSTANTE PAR LA METHODE *
! DES MOINDRES CARRES D UNE DROITE Y = A * X + B          *
!**********************************************************
!
! IN  X     : COORDONNEES DES POINTS (RAYON)
! IN  Y     : COORDONNEES DES POINTS (ERREUR MOYENNE)
! IN  NBPTS : NOMBRE DE POINTS
! OUT PE    : DEGRE DE LA SINGULARITE
!
    implicit none
!
! DECLARATION GLOBALE
!
    integer(kind=8) :: nbpts
    real(kind=8) :: x(nbpts), y(nbpts), pe
!
! DECLARATION LOCALE
!
    integer(kind=8) :: ipt, iter
    real(kind=8) :: xmoy, ymoy, vx, vxy, a, b, c, c1, f, ymin
!
! 1 - INITIALISATION
!
    c = 0.d+0
    iter = 0
!
    ymin = y(1)
    do ipt = 2, nbpts
        ymin = min(ymin, y(ipt))
    end do
!
100 continue
!
    xmoy = 0.d0
    ymoy = 0.d0
    vx = 0.d0
    vxy = 0.d0
!
! 2 - RECHERCHE DU X MOYEN ET DU Y MOYEN
!
    do ipt = 1, nbpts
        xmoy = xmoy+log(x(ipt))
        ymoy = ymoy+log(y(ipt)-c)
    end do
!
    xmoy = xmoy/nbpts
    ymoy = ymoy/nbpts
!
    do ipt = 1, nbpts
        vx = vx+(xmoy-log(x(ipt)))*(xmoy-log(x(ipt)))
        vxy = vxy+(xmoy-log(x(ipt)))*(ymoy-log(y(ipt)-c))
    end do
!
! 3 - PENTE THEORIQUE = VXY/VX
! CONSTANTE = Y - PENTE * X
!
    a = vxy/vx
    b = exp(ymoy-a*xmoy)
    f = 0.d+0
!
    do ipt = 1, nbpts
        f = f+(y(ipt)-c-b*x(ipt)**a)**2
    end do
!
    c1 = 0.d+0
!
    do ipt = 1, nbpts
        c1 = c1+(y(ipt)-c-b*x(ipt)**a)
    end do
!
    c1 = c1/nbpts
!
    if (c1 .gt. 0 .and. c1 .lt. ymin) then
        if (abs((c1-c)/c1) .gt. 1) then
            c = c1
            f = 0.d+0
            do ipt = 1, nbpts
                f = f+(y(ipt)-c-b*x(ipt)**a)**2
            end do
            iter = iter+1
            if (iter .lt. 5) goto 100
        end if
    end if
!
    pe = (a/2.d0)+1.d0
    if (pe .le. 0.4d0) pe = 0.4d0
!
end subroutine
