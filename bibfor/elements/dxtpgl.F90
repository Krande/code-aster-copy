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

subroutine dxtpgl(xyzg, pgl)
    implicit none
    real(kind=8) :: xyzg(3, *), pgl(3, 3)
!
!     IN  XYZG  R  9   COORDONNEES  X1 Y1 Z1 X2 Y2 ...
!     OUT PGL   R 3,3  MATRICE DE PASSAGE GLOBAL INTRINSEQUE
!     -----------------------------------------------------------------
!     CONSTRUCTION DE LA MATRICE DE PASSAGE GLOBAL --> INTRINSEQUE
!     POUR UNE MAILLE TRIANGLE DKT OU DST
!                                                        3
!                                                        *
!        I : VECTEUR UNITAIRE PORTE PAR 12              *  *
!                                                      *     *
!        K : PERPENDICULAIRE A 12 ET A 13             *        *
!                                                    *           *
!        J : PRODUIT VECTORIEL K I                  ****************
!                                                  1               2
!
!     -----------------------------------------------------------------
    real(kind=8) :: x21, y21, z21, x31, y31, z31, norm
!     -----------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    x21 = xyzg(1, 2)-xyzg(1, 1)
    y21 = xyzg(2, 2)-xyzg(2, 1)
    z21 = xyzg(3, 2)-xyzg(3, 1)
    x31 = xyzg(1, 3)-xyzg(1, 1)
    y31 = xyzg(2, 3)-xyzg(2, 1)
    z31 = xyzg(3, 3)-xyzg(3, 1)
!
    norm = sqrt(x21*x21+y21*y21+z21*z21)
    pgl(1, 1) = x21/norm
    pgl(1, 2) = y21/norm
    pgl(1, 3) = z21/norm
!
    pgl(3, 1) = y21*z31-z21*y31
    pgl(3, 2) = z21*x31-x21*z31
    pgl(3, 3) = x21*y31-y21*x31
!
    norm = sqrt(pgl(3, 1)*pgl(3, 1)+pgl(3, 2)*pgl(3, 2)+pgl(3, 3)*pgl(3, 3))
    pgl(3, 1) = pgl(3, 1)/norm
    pgl(3, 2) = pgl(3, 2)/norm
    pgl(3, 3) = pgl(3, 3)/norm
!
    pgl(2, 1) = pgl(3, 2)*pgl(1, 3)-pgl(3, 3)*pgl(1, 2)
    pgl(2, 2) = pgl(3, 3)*pgl(1, 1)-pgl(3, 1)*pgl(1, 3)
    pgl(2, 3) = pgl(3, 1)*pgl(1, 2)-pgl(3, 2)*pgl(1, 1)
!
end subroutine
