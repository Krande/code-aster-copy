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
subroutine gdmb(ne, kp, ajacob, en, enprim, &
                x0pg, b)
!
! FONCTION: POUR UN ELEMENT DE POUTRE EN GRAND DEPLACEMENT, CALCULE LA
!           CONTRIBUTION DU DEPLACEMENT DU NOEUD NE A LA MATRICE DE
!           DEFORMATION B AU POINT DE GAUSS KP.
!
!     IN  : NE        : NUMERO DU NOEUD
!           KP        : NUMERO DU POINT DE GAUSS
!           AJACOB    : JACOBIEN
!           EN        : FONCTIONS DE FORME
!           ENPRIM    : DERIVEES DES FONCTIONS DE FORME
!           X0PG      : DERIVEES DES COORDONNEES PAR RAP. A L'ABS. CURV.
!
!     OUT : B         : MATRICE DE DEFORMATION 6*6
! ------------------------------------------------------------------
    implicit none
#include "asterfort/antisy.h"
    real(kind=8) :: en(3, 2), enprim(3, 2), x0pg(3), b(6, 6), amat(3, 3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: kp, l, m, ne
    real(kind=8) :: ajacob, form, formpr, un, unsurj, zero
!-----------------------------------------------------------------------
    zero = 0.d0
    un = 1.d0
    do m = 1, 6
        do l = 1, 6
            b(l, m) = zero
        end do
    end do
    unsurj = un/ajacob
    form = en(ne, kp)
    formpr = unsurj*enprim(ne, kp)
    do l = 1, 6
        b(l, l) = formpr
    end do
    call antisy(x0pg, un, amat)
    do m = 1, 3
        do l = 1, 3
            b(l, m+3) = form*amat(l, m)
        end do
    end do
end subroutine
