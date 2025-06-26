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
subroutine abscvf(ndim, tabar, xe, s)
! aslint: disable=W1306
    implicit none
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/fcthyp.h"
    real(kind=8) :: xe, s, tabar(*)
    integer(kind=8) :: ndim
!
!
!                      TROUVER L'ABSCISSE CURVILIGNE D'UN POINT
!                      SUR UNE ARETE QUADRATIQUE A PARTIR DE SES
!                      COORDONNEES DANS L'ELEMENT DE REFERENCE
!
!     ENTREE
!       NDIM    : DIMENSION TOPOLOGIQUE DU MAILLAGE
!       TABAR  : COORDONNEES DES 3 NOEUDS QUI DEFINISSENT L'ARETE
!       XE     : COORDONNEES DU POINT DANS L'ELEMENT DE REFERENCE
!
!     SORTIE
!       S        : ABSCISSE CURVILIGNE DU POINT SUR L'ARETE
!......................................................................
!
    real(kind=8) :: coef1, coef2, coef3, coef4, eps
    real(kind=8) :: pt1(ndim), pt2(ndim), pt3(ndim)
    real(kind=8) :: d, mu
    real(kind=8) :: borsup, borinf
    integer(kind=8) :: i
    character(len=8) :: typfct
!
!......................................................................
!
!
!     TABAR : XE2=-1  /  XE1= 1  /  XE3= 0
!     XE2 XENT LE POINT D'ORIGINE
!
!     CALCUL DE COEF1, COEF2, COEF3, D
    coef1 = 0.d0
    coef2 = 0.d0
    coef3 = 0.d0
    pt1(:) = 0.d0
    pt2(:) = 0.d0
    pt3(:) = 0.d0
    eps = 1.d-7
!
    do i = 1, ndim
        pt1(i) = tabar(i)
        pt2(i) = tabar(ndim+i)
        pt3(i) = tabar(2*ndim+i)
    end do
!
    do i = 1, ndim
        coef1 = coef1+(pt1(i)-2*pt3(i)+pt2(i))*(pt1(i)-2*pt3(i)+pt2(i))
    end do
!
    do i = 1, ndim
        coef2 = coef2+(pt2(i)-pt1(i))*(pt1(i)-2*pt3(i)+pt2(i))
    end do
!
    do i = 1, ndim
        coef3 = coef3+(pt2(i)-pt1(i))*(pt2(i)-pt1(i))/4
    end do
!
    d = coef2*coef2-4*coef1*coef3
!
!     CALCUL ABSCISSE CURVILIGNE DU POINT
!
    if (abs(coef1) .le. eps) then
        s = (xe+1)*sqrt(coef3)
    else if (abs(coef1) .gt. r8prem()) then
        if (abs(d) .le. r8prem()) then
            s = (coef1*xe*xe+coef2*xe+coef2-coef1)/(2*sqrt(coef1))
        else if (d .gt. eps) then
            mu = sqrt(d/(4*coef1*coef1))
            coef4 = mu*mu*sqrt(coef1)/4
            typfct = 'ACOSH'
            call fcthyp(typfct, (2*coef1*xe+coef2)/(2*coef1*mu), borsup)
            call fcthyp(typfct, (coef2-2*coef1)/(2*coef1*mu), borinf)
            s = coef4*(sinh(2*borsup)-2*borsup)-coef4*(sinh(2*borinf)-2*borinf)
        else if (d .lt. eps) then
            mu = sqrt(-d/(4*coef1*coef1))
            coef4 = mu*mu*sqrt(coef1)/4
            typfct = 'ASINH'
            call fcthyp(typfct, (2*coef1*xe+coef2)/(2*coef1*mu), borsup)
            call fcthyp(typfct, (coef2-2*coef1)/(2*coef1*mu), borinf)
            s = coef4*(sinh(2*borsup)+2*borsup)-coef4*(sinh(2*borinf)+2*borinf)
!
        end if
    end if
!
    if (s .lt. 0.d0) then
        ASSERT(2 .eq. 3)
    end if
!
end subroutine
