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
subroutine gdmmas(kp, nno, pjacob, en, grani, &
                  rot0, mass)
!
! FONCTION: POUR UN ELEMENT DE POUTRE EN GRAND DEPLACEMENT, CALCULE LA
!           MATRICE DE MASSE EN POSITION DE REFERENCE.
!
!     IN  : KP        : NUMERO DU POINT DE GAUSS
!           NNO       : NOMBRE DE NOEUDS
!           PJACOB    : POIDS * JACOBIEN
!           EN        : FONCTIONS DE FORME
!           GRANI     : DIAGONALE DU TENSEUR D'INERTIE EN AXES LOCAUX
!                       POUR LES 3 1ERES COMPOSANTES, RHO*A POUR LA 4EME
!           ROT0      : MATRICE DE ROTATION DES AXES PRINCIPAUX D'INERT.
!                       AU POINT DE GAUSS DANS LA POSITION DE REFERENCE,
!                       PAR RAPPORT AUX AXES GENERAUX
!
!     OUT : MASS      : MATRICE DE MASSE (CUMUL DES CONTRIBUTIONS DES
!                       POINTS DE GAUSS)
! ------------------------------------------------------------------
    implicit none
#include "asterfort/cumuma.h"
#include "asterfort/promat.h"
    real(kind=8) :: en(3, 2), grani(4), rot0(3, 3), mass(18, 18), imas(6, 6)
    real(kind=8) :: iro(3, 3), amat1(3, 3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, kp, nno
    real(kind=8) :: coef, pjacob, zero
!-----------------------------------------------------------------------
    zero = 0.d0
    do j = 1, 6
        do i = 1, 6
            imas(i, j) = zero
        end do
    end do
!
    do i = 1, 3
        imas(i, i) = grani(4)
    end do
!
    do j = 1, 3
        do i = 1, 3
            amat1(i, j) = grani(i)*rot0(j, i)
        end do
    end do
    call promat(rot0, 3, 3, 3, amat1, &
                3, 3, 3, iro)
!
    do j = 1, 3
        do i = 1, 3
            imas(3+i, 3+j) = iro(i, j)
        end do
    end do
!
    do j = 1, nno
        do i = 1, nno
            coef = pjacob*en(i, kp)*en(j, kp)
            call cumuma(i, j, imas, coef, mass)
        end do
    end do
end subroutine
