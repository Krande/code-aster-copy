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

subroutine calc_h_tria(ino, x3d1, x3d2, x3d3, h)
    implicit none
!  DESCRIPTION :
!  -----------
!       CALCUL DE LA HAUTEUR H RELATIVE AU NOEUD INO DU TRIANGLE DEFINI
!       PAR LES POINTS X3D1, X3D2, X3D3
!
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
#include "asterfort/normev.h"
! ARGUMENTS
! ---------
    integer(kind=8) :: ino
    real(kind=8) :: x3d1(3), x3d2(3), x3d3(3), h
!
! VARIABLES LOCALES
! -----------------
    real(kind=8) :: x3d(3, 3), vect_jk(3), norme_jk, proj, vect_ji(3)
    integer(kind=8) :: jno, kno
!
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
!
    x3d(1, 1) = x3d1(1)
    x3d(2, 1) = x3d1(2)
    x3d(3, 1) = x3d1(3)
!
    x3d(1, 2) = x3d2(1)
    x3d(2, 2) = x3d2(2)
    x3d(3, 2) = x3d2(3)
!
    x3d(1, 3) = x3d3(1)
    x3d(2, 3) = x3d3(2)
    x3d(3, 3) = x3d3(3)
!
    jno = mod(ino+1, 3)
    if (jno .eq. 0) jno = 3
    kno = mod(ino+2, 3)
    if (kno .eq. 0) kno = 3
!
!   vecteur directeur du coté opposé au noeud
    vect_jk(1) = x3d(1, kno)-x3d(1, jno)
    vect_jk(2) = x3d(2, kno)-x3d(2, jno)
    vect_jk(3) = x3d(3, kno)-x3d(3, jno)
!
    call normev(vect_jk, norme_jk)
!
    vect_ji(1) = x3d(1, ino)-x3d(1, jno)
    vect_ji(2) = x3d(2, ino)-x3d(2, jno)
    vect_ji(3) = x3d(3, ino)-x3d(3, jno)
!
    proj = (vect_jk(1)*vect_ji(1)+vect_jk(2)*vect_ji(2)+vect_jk(3)*vect_ji(3))
!
    vect_ji(1) = vect_ji(1)-proj*vect_jk(1)
    vect_ji(2) = vect_ji(2)-proj*vect_jk(2)
    vect_ji(3) = vect_ji(3)-proj*vect_jk(3)
!
    call normev(vect_ji, h)

end subroutine
