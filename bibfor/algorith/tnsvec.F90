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
subroutine tnsvec(choix, ndim, mat, vec, r)
    implicit none
#include "asterfort/assert.h"
    integer(kind=8) :: ndim, i, choix
    real(kind=8) :: mat(3, 3), vec(2*ndim), r
!     FONCTION TRANSFORMANT UNE MATRICE SYMETRIQUE (TENSEUR) 3,3
!     EN VECTEUR 6 (OU 4 EN 2D) ET MULTIPLIANT LES TERMES NON DIAGONAUX
!     PAR UN REEL R ET INVERSEMENT.
!     IN : CHOIX =3 OU 6 (DIMENSION DE L'OBJET EN ENTREE)
!     Si CHOIX=3
!        IN : MAT
!        OUT : VEC
!     SI CHOIX=6
!        IN : VEC
!        OUT : MAT
!
    if (choix .eq. 3) then
!
!        TRANSFORMATION MATRICE EN VECTEUR
        do i = 1, 3
            vec(i) = mat(i, i)
        end do
        vec(4) = mat(1, 2)*r
        if (ndim .eq. 3) then
            vec(5) = mat(1, 3)*r
            vec(6) = mat(2, 3)*r
        end if
!
!
    else if (choix .eq. 6) then
!
!        TRANSFORMATION VECTEUR EN MATRICE
        do i = 1, 3
            mat(i, i) = vec(i)
        end do
        mat(1, 2) = vec(4)*r
        mat(2, 1) = vec(4)*r
        if (ndim .eq. 2) then
            mat(1, 3) = 0.d0
            mat(2, 3) = 0.d0
            mat(3, 1) = 0.d0
            mat(3, 2) = 0.d0
        else
            mat(1, 3) = vec(5)*r
            mat(3, 1) = vec(5)*r
            mat(2, 3) = vec(6)*r
            mat(3, 2) = vec(6)*r
        end if
    else
        ASSERT(choix .eq. 3)
    end if
!
end subroutine
