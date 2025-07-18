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

subroutine pj2da2(ino2, geom2, i, geom1, tria3, &
                  cobary, d2, surf)
    implicit none
#include "asterc/r8maem.h"
    real(kind=8) :: cobary(3), geom1(*), geom2(*), d2, surf
    integer(kind=8) :: ino2, i, tria3(*)

!  but :
!    determiner la distance d2 entre le noeud ino2 et le tria3 i.
!    determiner les coordonnees barycentriques
!    du point de i le plus proche de ino2.

!  in   ino2       i  : numero du noeud de m2 cherche
!  in   geom2(*)   r  : coordonnees des noeuds du maillage m2
!  in   geom1(*)   r  : coordonnees des noeuds du maillage m1
!  in   i          i  : numero du tria3 candidat
!  in   tria3(*)   i  : objet '&&pjxxco.tria3'
!  out  cobary(3)  r  : coordonnees barycentriques de ino2 projete sur i
!  out  d2         r  : carre de la distance entre i et ino2
!  out  surf       r  : surface du tria3 i

! ----------------------------------------------------------------------
    real(kind=8) :: x1, y1, x2, y2, x3, y3, xp, yp, xm, ym
    real(kind=8) :: ksi, dist
    real(kind=8) :: v1(2), v2(2), v3(2), m(2), xc
! DEB ------------------------------------------------------------------
    xm = geom2(3*(ino2-1)+1)
    ym = geom2(3*(ino2-1)+2)

    x1 = geom1(3*(tria3(1+4*(i-1)+1)-1)+1)
    y1 = geom1(3*(tria3(1+4*(i-1)+1)-1)+2)
    x2 = geom1(3*(tria3(1+4*(i-1)+2)-1)+1)
    y2 = geom1(3*(tria3(1+4*(i-1)+2)-1)+2)
    x3 = geom1(3*(tria3(1+4*(i-1)+3)-1)+1)
    y3 = geom1(3*(tria3(1+4*(i-1)+3)-1)+2)

    v1(1) = x3-x2
    v1(2) = y3-y2
    v2(1) = x1-x3
    v2(2) = y1-y3
    v3(1) = x2-x1
    v3(2) = y2-y1

    surf = abs(v2(1)*v3(2)-v2(2)*v3(1))
    surf = surf/2.d0
    d2 = r8maem()

!     COTE 1 (2->3):
!     --------------
    m(1) = xm-x2
    m(2) = ym-y2
    xc = v1(1)*v1(1)+v1(2)*v1(2)
    if (xc .eq. 0) then
        ksi = 0.5d0
    else
        ksi = (m(1)*v1(1)+m(2)*v1(2))/xc
    end if
    if (ksi .ge. 1.d0) then
        ksi = 1.d0
    else if (ksi .le. 0.d0) then
        ksi = 0.d0
    end if
    xp = ksi*x3+(1.d0-ksi)*x2
    yp = ksi*y3+(1.d0-ksi)*y2
    dist = (xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)
    if (dist .lt. d2) then
        cobary(1) = 0.d0
        cobary(2) = 1.d0-ksi
        cobary(3) = ksi
        d2 = dist
    end if

!     COTE 2 (3->1):
!     --------------
    m(1) = xm-x3
    m(2) = ym-y3
    xc = v2(1)*v2(1)+v2(2)*v2(2)
    if (xc .eq. 0) then
        ksi = 0.5d0
    else
        ksi = (m(1)*v2(1)+m(2)*v2(2))/xc
    end if
    if (ksi .ge. 1.d0) then
        ksi = 1.d0
    else if (ksi .le. 0.d0) then
        ksi = 0.d0
    end if
    xp = ksi*x1+(1.d0-ksi)*x3
    yp = ksi*y1+(1.d0-ksi)*y3
    dist = (xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)
    if (dist .lt. d2) then
        cobary(2) = 0.d0
        cobary(3) = 1.d0-ksi
        cobary(1) = ksi
        d2 = dist
    end if

!     COTE 3 (1->2):
!     --------------
    m(1) = xm-x1
    m(2) = ym-y1
    xc = v3(1)*v3(1)+v3(2)*v3(2)
    if (xc .eq. 0) then
        ksi = 0.5d0
    else
        ksi = (m(1)*v3(1)+m(2)*v3(2))/xc
    end if
    if (ksi .ge. 1.d0) then
        ksi = 1.d0
    else if (ksi .le. 0.d0) then
        ksi = 0.d0
    end if
    xp = ksi*x2+(1.d0-ksi)*x1
    yp = ksi*y2+(1.d0-ksi)*y1
    dist = (xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)
    if (dist .lt. d2) then
        cobary(3) = 0.d0
        cobary(1) = 1.d0-ksi
        cobary(2) = ksi
        d2 = dist
    end if

end subroutine
