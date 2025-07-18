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
subroutine pj3dgb(ino2, geom2, geom1, tetr4, ndec, &
                  btdi, btvr, btnb, btlc, btco, &
                  p1, q1, r1, p2, q2, &
                  r2)
    implicit none
#include "asterfort/assert.h"
    real(kind=8) :: geom1(*), geom2(*), btvr(*)
    integer(kind=8) :: ino2, p1, q1, r1, p2, q2, r2, ndec
    integer(kind=8) :: btdi(*), btnb(*), btlc(*), btco(*), tetr4(*)
!     BUT :
!       TROUVER LA "GROSSE BOITE" (P1,Q1,R1,P2,Q2,R2) DANS LA QUELLE
!       ON EST SUR DE TROUVER LE TETR4/TRIA3/SEG2 LE PLUS PROCHE DE INO2
!
!  IN   INO2       I  : NUMERO DU NOEUD DE M2 CHERCHE
!  IN   GEOM2(*)   R  : COORDONNEES DES NOEUDS DU MAILLAGE M2
!  IN   GEOM1(*)   R  : COORDONNEES DES NOEUDS DU MAILLAGE M1
!  IN   TETR4(*)   I  : OBJET '&&PJXXCO.TETR4/TRIA3/SEG2'
!  IN   NDEC       I  : 6  SI '&&PJXXCO.TETR4'
!                     : 4  SI '&&PJXXCO.TRIA3'
!                     : 3  SI '&&PJXXCO.SEG2'
!  IN   BTDI(*)    I  : OBJET .BT3DDI DE LA SD BOITE_3D
!  IN   BTVR(*)    R  : OBJET .BT3DVR DE LA SD BOITE_3D
!  IN   BTNB(*)    I  : OBJET .BT3DNB DE LA SD BOITE_3D
!  IN   BTLC(*)    I  : OBJET .BT3DLC DE LA SD BOITE_3D
!  IN   BTCO(*)    I  : OBJET .BT3DCO DE LA SD BOITE_3D
!  OUT  P1         I  : ABSCISSE DU COIN BAS/GAUCHE DE LA GROSSE BOITE
!  OUT  Q1         I  : ORDONNEE DU COIN BAS/GAUCHE DE LA GROSSE BOITE
!  OUT  R1         I  : COTE     DU COIN BAS/GAUCHE DE LA GROSSE BOITE
!  OUT  P2         I  : ABSCISSE DU COIN HAUT/DROIT DE LA GROSSE BOITE
!  OUT  Q2         I  : ORDONNEE DU COIN HAUT/DROIT DE LA GROSSE BOITE
!  OUT  R2         I  : COTE     DU COIN HAUT/DROIT DE LA GROSSE BOITE
! ----------------------------------------------------------------------
    real(kind=8) :: d, x1, y1, z1, x2, y2, z2, xmin, ymin, zmin, dx, dy, dz
    integer(kind=8) :: p0, q0, r0, nx, ny, nz, k, p, q, r, itr, ntrbt, iposi, ino1
!
! DEB ------------------------------------------------------------------
!
    nx = btdi(1)
    ny = btdi(2)
    nz = btdi(3)
    dx = btvr(7)
    dy = btvr(8)
    dz = btvr(9)
    xmin = btvr(1)
    ymin = btvr(3)
    zmin = btvr(5)
!
!
!      1. ON CHERCHE UNE BOITE NON VIDE AUTOUR DE INO2
!     -------------------------------------------------------
    p0 = int((geom2(3*(ino2-1)+1)-xmin)/dx)+1
    q0 = int((geom2(3*(ino2-1)+2)-ymin)/dy)+1
    r0 = int((geom2(3*(ino2-1)+3)-zmin)/dz)+1
!
    do k = 0, max(nx, ny, nz)-1
        do p = max(p0-k, 1), min(p0+k, nx)
            do q = max(q0-k, 1), min(q0+k, ny)
                do r = max(r0-k, 1), min(r0+k, nz)
                    ntrbt = btnb((r-1)*nx*ny+(q-1)*nx+p)
!             -- SI LA BOITE EST NON VIDE :
                    if (ntrbt .gt. 0) then
!               -- ON CHOISIT LE 1ER NOEUD DU 1ER TETR4 DE LA BOITE:INO1
                        iposi = btlc((r-1)*nx*ny+(q-1)*nx+p)
                        itr = btco(iposi+1)
                        ino1 = tetr4(1+ndec*(itr-1)+1)
                        goto 50
                    end if
                end do
            end do
        end do
    end do
    ASSERT(.false.)
!
!
50  continue
!     2. ON CALCULE LA DISTANCE ENTRE INO2 ET INO1
!     -------------------------------------------------------
    x1 = geom1(3*(ino1-1)+1)
    y1 = geom1(3*(ino1-1)+2)
    z1 = geom1(3*(ino1-1)+3)
    x2 = geom2(3*(ino2-1)+1)
    y2 = geom2(3*(ino2-1)+2)
    z2 = geom2(3*(ino2-1)+3)
    d = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
!
!
!     3. ON DETERMINE LA GROSSE BOITE CONTENANT :
!        INO2 - D*VECTEUR_I - D*VECTEUR_J - D*VECTEUR_K
!     ET INO2 + D*VECTEUR_I + D*VECTEUR_J + D*VECTEUR_K
!     -------------------------------------------------------
    p1 = int((x2-d-xmin)/dx)+1
    q1 = int((y2-d-ymin)/dy)+1
    r1 = int((z2-d-zmin)/dz)+1
    p1 = max(1, p1)
    q1 = max(1, q1)
    r1 = max(1, r1)
!
    p2 = int((x2+d-xmin)/dx)+1
    q2 = int((y2+d-ymin)/dy)+1
    r2 = int((z2+d-zmin)/dz)+1
    p2 = min(nx, p2)
    q2 = min(ny, q2)
    r2 = min(nz, r2)
!
!
end subroutine
