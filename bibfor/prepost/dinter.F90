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

subroutine dinter(coorc, ray, coor1, coor2, coorin)
!
!*********************************************************
!              BUT DE CETTE ROUTINE :                    *
! CALCULER LES COORDONNEES DU POINT D INTERSECTION ENTRE *
! UNE DROITE ET UN CERCLE DE CENTRE X ET DE RAYON RAY    *
!*********************************************************
!
! IN  COORC(2)   : COORDONNEES DU CENTRE X
! IN  RAY        : RAYON DU CERCLE
! IN  COOR1(2)   : COORDONNEES DU NOEUD 1 DE LA DROITE
! IN  COOR2(2)   : COORDONNEES DU NOEUD 2 DE LA DROITE
! OUT COORIN(2) : COORDONNEES DE L INTERSECTION
!
    implicit none
#include "asterfort/utmess.h"
!
! DECLARATION GLOBALE
!
    real(kind=8) :: coorc(2), ray, coor1(2), coor2(2), coorin(2)
!
! DECLARATION LOCALE
!
    real(kind=8) :: xc, yc, x1, y1, x2, y2, xi, yi
    real(kind=8) :: det, val, m
    real(kind=8) :: prec
    parameter(prec=10.d-6)
!
    xc = coorc(1)
    yc = coorc(2)
    x1 = coor1(1)
    y1 = coor1(2)
    x2 = coor2(1)
    y2 = coor2(2)
!
    if ((y2-y1) .ge. prec .or. (y2-y1) .le. -prec) then
        m = (x2-x1)/(y2-y1)
        val = x1-xc-m*y1
        det = (yc-m*val)**2-(1.d0+m**2)*(val**2+yc**2-ray**2)
        if (det .ge. (-prec)) then
            det = abs(det)
            yi = (yc-m*val+sqrt(det))/(1.d0+m**2)
            if (((yi-y1) .le. prec .and. (yi-y2) .ge. -prec) .or. &
                ((yi-y1) .ge. -prec .and. (yi-y2) .le. prec)) then
                coorin(2) = yi
                coorin(1) = x1+m*(yi-y1)
            else
                yi = (yc-m*val-sqrt(det))/(1.d0+m**2)
                if (((yi-y1) .le. prec .and. (yi-y2) .ge. -prec) .or. &
                    ((yi-y1) .ge. -prec .and. (yi-y2) .le. prec)) then
                    coorin(2) = yi
                    coorin(1) = x1+m*(yi-y1)
                else
                    call utmess('F', 'PREPOST_21')
                end if
            end if
        else
!          R1 = SQRT((XC-X1)**2 +(YC-Y1)**2 )
!          R2 = SQRT((XC-X2)**2 +(YC-Y2)**2 )
        end if
    else
        det = ray**2-(y1-yc)**2
        if (det .ge. (-prec)) then
            det = abs(det)
            xi = xc+sqrt(det)
            if (((xi-x1) .le. prec .and. (xi-x2) .ge. -prec) .or. &
                ((xi-x1) .ge. -prec .and. (xi-x2) .le. prec)) then
                coorin(1) = xi
                coorin(2) = y1
            else
                xi = xc-sqrt(det)
                if (((xi-x1) .le. prec .and. (xi-x2) .ge. -prec) .or. &
                    ((xi-x1) .ge. -prec .and. (xi-x2) .le. prec)) then
                    coorin(1) = xi
                    coorin(2) = y1
                else
                    call utmess('F', 'PREPOST_21')
                end if
            end if
        else
            call utmess('F', 'PREPOST_21')
        end if
    end if
end subroutine
