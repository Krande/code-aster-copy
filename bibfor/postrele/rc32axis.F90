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

subroutine rc32axis(nbabsc, absc, xcoo, ycoo, vale, momen0_axis, momen1_axis, momen2_axis, rho)
    implicit none
#include "asterfort/rcrot.h"
#include "asterc/r8pi.h"
    integer(kind=8) :: nbabsc
    real(kind=8) :: absc(nbabsc), vale(4, nbabsc), xcoo(nbabsc)
    real(kind=8) :: ycoo(nbabsc), momen0_axis(4), momen1_axis(4), momen2_axis(4), rho
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3200
!
!     METHODE PRSECT DE ANSYS POUR LINEARISER LES CONTRAINTES
!     PRISE EN COMPTE DE LA DISTRIBUTION RADIALE INHOMOGENE DU MATERIAU D'UN MODELE AXISYMETRIQUE
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbsgt
    real(kind=8) :: sigm(4, nbabsc), f(nbabsc)
    real(kind=8) :: l, x1, x2, y1, y2, phi, pi, mid_rad, xf, xh, dx
! DEB ------------------------------------------------------------------
!
!
    pi = r8pi()
    l = absc(nbabsc)-absc(1)
    nbsgt = nbabsc-1
    dx = l/nbsgt
    x1 = xcoo(1)
    y1 = ycoo(1)
    x2 = xcoo(nbabsc)
    y2 = ycoo(nbabsc)
    if (x1 .eq. x2) then
        phi = pi/2
    else
        phi = atan((y1-y2)/(x1-x2))
    end if
    call rcrot(nbabsc, phi, vale, sigm)
    mid_rad = (x1+x2)/2
    xf = l**2*cos(phi)/(12*mid_rad)
!
! --- CALCULER LES CONTRAINTES DE MEMBRANE
! --- sigm_x^m
    momen0_axis(1) = (sigm(1, 1)/2+sigm(1, nbabsc)/2+sum(sigm(1, 2:nbabsc-1)))/nbsgt
! --- sigm_y^m
    f = sigm(2, :)*xcoo
    momen0_axis(2) = (f(1)/2+f(nbabsc)/2+sum(f(2:nbabsc-1)))*dx/(mid_rad*l)
! --- sigm_z^m
    if (rho .gt. 0.d0) then
        f = sigm(3, :)*(1+xcoo/rho)
        momen0_axis(3) = (f(1)/2+f(nbabsc)/2+sum(f(2:nbabsc-1)))/nbsgt
    else
        momen0_axis(3) = (sigm(3, 1)/2+sigm(3, nbabsc)/2+sum(sigm(3, 2:nbabsc-1)))/nbsgt
    end if
! --- sigm_xy^m
    f = sigm(4, :)*xcoo
    momen0_axis(4) = (f(1)/2+f(nbabsc)/2+sum(f(2:nbabsc-1)))*dx/(mid_rad*l)
!
! --- CALCULER LES CONTRAINTES DE FLEXION
! --- sigm_x^b (kbr = 0)
    momen1_axis(1) = sigm(1, 1)-momen0_axis(1)
    momen2_axis(1) = sigm(1, nbabsc)-momen0_axis(1)
! --- sigm_y^b
    f = (absc-l/2-xf)*sigm(2, :)*xcoo
    momen1_axis(2) = (-l/2-xf)/(mid_rad*l*(l**2/12-xf**2))*(f(1)/2+ &
                                                            f(nbabsc)/2+sum(f(2:nbabsc-1)))*dx
    momen2_axis(2) = momen1_axis(2)*(xf-l/2)/(xf+l/2)
! --- sigm_z^b
    if (rho .gt. 0.d0) then
        xh = l**2/(12*rho)
        f = (absc-l/2-xh)*sigm(3, :)*(1+(absc-l/2)/rho)
        momen1_axis(3) = (-l/2-xh)/(l*(l**2/12-xh**2))*(f(1)/2+f(nbabsc)/2+sum(f(2:nbabsc-1)))*dx
        momen2_axis(3) = momen1_axis(3)*((xh-l/2)/(xh+l/2))
    else
        f = (absc-l/2)*sigm(3, :)
        momen1_axis(3) = -6/l**2*(f(1)/2+f(nbabsc)/2+sum(f(2:nbabsc-1)))*dx
        momen2_axis(3) = -momen1_axis(3)
    end if
! --- sigm_xy^b = 0
    momen1_axis(4) = 0.d0
    momen2_axis(4) = 0.d0
end subroutine
!
