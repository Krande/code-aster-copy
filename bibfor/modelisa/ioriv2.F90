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
function ioriv2(num, n, noeud, vect, coor)
!.======================================================================
    implicit none
#include "asterfort/utmess.h"
!
!     IORIV2  --  ORIENTATION D'UNE MAILLE PAR RAPPORT A UN VECTEUR
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NUM          IN/OUT  K*     NUMEROTATION DE LA MAILLE
!    N              IN    K*     NOMBRE DE NOEUDS DE LA MAILLE
!
!   CODE RETOUR IORIV2 : 0 SI LA MAILLE NE CONTIENT PAS LE NOEUD
!                       -1 OU 1 SINON (SELON QU'IL AIT OU NON
!                                      FALLU REORIENTER)
    integer(kind=8) :: n, num(n)
    real(kind=8) :: vect(3), coor(3, *)
!
!     DONNEES POUR TRIA3,TRIA6,TRIA7,QUAD4,QUAD8,QUAD9
!     NOMBRE DE SOMMETS EN FONCTION DU NOMBRE DE NOEUDS DE L'ELEMENT
    integer(kind=8) :: nsom(9)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i2, i3, ioriv2, k, l
    integer(kind=8) :: n1, n2, n3, noeud, nso
    real(kind=8) :: scal, x1, x2, x3, xa, xb, xn
    real(kind=8) :: y1, y2, y3, ya, yb, yn, z1
    real(kind=8) :: z2, z3, za, zb, zn
!-----------------------------------------------------------------------
    data nsom/0, 0, 3, 4, 0, 3, 3, 4, 4/
!
#define x(i) coor(1,i)
#define y(i) coor(2,i)
#define z(i) coor(3,i)
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    nso = nsom(n)
!     BOUCLE SUR LES SOMMETS
    ioriv2 = 0
    do i = 1, nso
        if (num(i) .eq. noeud) then
            i1 = i
            i2 = i-1
            if (i .eq. 1) i2 = nso
            i3 = i+1
            if (i .eq. nso) i3 = 1
            n1 = num(i1)
            n2 = num(i2)
            n3 = num(i3)
            x1 = x(n1)
            y1 = y(n1)
            z1 = z(n1)
            x2 = x(n2)
            y2 = y(n2)
            z2 = z(n2)
            x3 = x(n3)
            y3 = y(n3)
            z3 = z(n3)
            xa = x2-x1
            ya = y2-y1
            za = z2-z1
            xb = x3-x1
            yb = y3-y1
            zb = z3-z1
!     VECTEUR NORMAL AU PLAN TANGENT AU NOEUD
            xn = ya*zb-yb*za
            yn = za*xb-zb*xa
            zn = xa*yb-xb*ya
            scal = xn*vect(1)+yn*vect(2)+zn*vect(3)
            if (scal .lt. 0) then
                ioriv2 = 1
            else if (scal .gt. 0) then
                ioriv2 = -1
            else
                call utmess('F', 'MODELISA4_76')
            end if
        end if
    end do
    if (ioriv2 .lt. 0) then
!       ON PERMUTE LES SOMMETS
        k = num(2)
        l = num(nso)
        num(2) = l
        num(nso) = k
!       ON PERMUTE LES NOEUDS INTERMEDIAIRES
        if (n .ne. nso) then
            do i = 1, nso/2
                k = num(nso+i)
                l = num(2*nso+1-i)
                num(nso+i) = l
                num(2*nso+1-i) = k
            end do
        end if
    end if
!
end function
