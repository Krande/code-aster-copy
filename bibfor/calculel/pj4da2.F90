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
subroutine pj4da2(ino2, geom2, i, geom1, tria3, &
                  cobary, d2, surf)
    implicit none
#include "asterf_types.h"
#include "asterc/r8maem.h"
#include "asterfort/pj3da3.h"
#include "asterfort/pj3da4.h"
    real(kind=8) :: cobary(3), geom1(*), geom2(*), d2, surf
    integer(kind=8) :: ino2, i, tria3(*)
!  but :
!    determiner la distance d entre le noeud ino2 et le tria3 i.
!    determiner les coordonnees barycentriques
!    du point de i le plus proche de ino2.
!
!  in   ino2       i  : numero du noeud de m2 cherche
!  in   geom2(*)   r  : coordonnees des noeuds du maillage m2
!  in   geom1(*)   r  : coordonnees des noeuds du maillage m1
!  in   i          i  : numero du tria3 candidat
!  in   tria3(*)   i  : objet '&&pjxxco.tria3'
!  out  cobary(3)  r  : coordonnees barycentriques de ino2 projete sur i
!  out  d2         r  : carre de la distance entre i et ino2
!  out  surf       r  : surface du tria3 i
!
! ----------------------------------------------------------------------
    integer(kind=8) :: k
    aster_logical :: ok
    real(kind=8) :: dp, l1, l2, l3, la, lb, lc
    real(kind=8) :: a(3), b(3), c(3), m(3), ab(3), ac(3), v(3)
! DEB ------------------------------------------------------------------
!
    do k = 1, 3
        m(k) = geom2(3*(ino2-1)+k)
        a(k) = geom1(3*(tria3(1+4*(i-1)+1)-1)+k)
        b(k) = geom1(3*(tria3(1+4*(i-1)+2)-1)+k)
        c(k) = geom1(3*(tria3(1+4*(i-1)+3)-1)+k)
        ab(k) = b(k)-a(k)
        ac(k) = c(k)-a(k)
    end do
!
    d2 = r8maem()
    dp = r8maem()
!
!
!   1. on cherche le point le plus proche a l'interieur du tria3
!   -------------------------------------------------------------
    call pj3da3(m, a, b, c, ok, &
                l1, l2, l3, dp)
    if ((ok) .and. (dp .lt. d2)) then
        d2 = dp
        la = l1
        lb = l2
        lc = l3
    end if
!
!
!   2. on boucle sur les 3 arretes du tria3 :
!   -----------------------------------------
    call pj3da4(m, a, b, l1, l2, &
                dp)
    if (dp .lt. d2) then
        d2 = dp
        la = l1
        lb = l2
        lc = 0.d0
    end if
!
    call pj3da4(m, b, c, l1, l2, &
                dp)
    if (dp .lt. d2) then
        d2 = dp
        lb = l1
        lc = l2
        la = 0.d0
    end if
!
    call pj3da4(m, a, c, l1, l2, &
                dp)
    if (dp .lt. d2) then
        d2 = dp
        la = l1
        lc = l2
        lb = 0.d0
    end if
!
!
!   3. on calcule surf :
!   --------------------
    v(1) = ab(2)*ac(3)-ab(3)*ac(2)
    v(2) = ab(3)*ac(1)-ab(1)*ac(3)
    v(3) = ab(1)*ac(2)-ab(2)*ac(1)
    surf = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
    surf = surf/2.d0
!
!
    cobary(1) = la
    cobary(2) = lb
    cobary(3) = lc
!
end subroutine
