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
subroutine gdfine(kp, nno, pjacob, en, grani, &
                  rot0, rotk, omgk, ompgk, fint)
!
! FONCTION: POUR UN ELEMENT DE POUTRE EN GRAND DEPLACEMENT, CALCULE LA
!           CONTRIBUTION DU POINT DE GAUSS NUMERO KP A LA PARTIE ROTA-
!           TOIRE DES FORCES D'INERTIE. CES FORCES D'INERTIE SONT RE-
!           TRANCHEES DES FORCES INTERNES.
!
!     IN  : KP        : NUMERO DU POINT DE GAUSS
!           NNO       : NOMBRE DE NOEUDS DE L'ELEMENT
!           PJACOB    : POIDS * JACOBIEN
!           EN        : FONCTIONS DE FORME
!           GRANI     : DIAGONALE DU TENSEUR D'INERTIE EN AXES LOCAUX
!                       POUR LES 3 1ERES COMPOSANTES, RHO*A POUR LA 4EME
!           ROT0      : MATRICE DE ROTATION DES AXES PRINCIPAUX D'INERT.
!                       AU POINT DE GAUSS DANS LA POSITION DE REFERENCE,
!                       PAR RAPPORT AUX AXES GENERAUX
!           ROTK      : MATRICE DE ROTATION ACTUELLE
!           OMGK      : VITESSE ANGULAIRE ACTUELLE
!           OMPGK     : ACCELERATION ANGULAIRE ACTUELLE
!
!     OUT : FINT      : FORCES INT. (CUMUL DES CONTRIB. DES PTS DE GAUS)
! ------------------------------------------------------------------
    implicit none
#include "asterfort/antisy.h"
#include "asterfort/promat.h"
    real(kind=8) :: irott(3, 3)
    real(kind=8) :: en(3, 2), grani(4), rot0(3, 3), rotk(3, 3), omgk(3)
    real(kind=8) :: ompgk(3), fint(6, 3), rotabs(3, 3), omegat(3, 3)
    real(kind=8) :: amati(3, 3), amat1(3, 3), amat2(3, 3), fors(6), v1(3), v2(3)
    real(kind=8) :: v3(3)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, kp, ne, nno
    real(kind=8) :: coef, pjacob, un
!-----------------------------------------------------------------------
    un = 1.d0
!
    call promat(rotk, 3, 3, 3, rot0, &
                3, 3, 3, rotabs)
!
    do j = 1, 3
        do i = 1, 3
            irott(i, j) = grani(i)*rotabs(j, i)
        end do
    end do
    call promat(rotabs, 3, 3, 3, irott, &
                3, 3, 3, amati)
    call promat(amati, 3, 3, 3, omgk, &
                3, 3, 1, v1)
!
    call antisy(omgk, un, omegat)
!
    call promat(omegat, 3, 3, 3, v1, &
                3, 3, 1, v2)
!
    call promat(amati, 3, 3, 3, ompgk, &
                3, 3, 1, v1)
!
!* ON CALCULE CE QU'APPORTE CHDYNL A LA PARTIE ROTATOIRE POUR LE
!* RETRANCHER ICI.
    do j = 1, 3
        do i = 1, 3
            amat1(i, j) = grani(i)*rot0(j, i)
        end do
    end do
    call promat(rot0, 3, 3, 3, amat1, &
                3, 3, 3, amat2)
    call promat(amat2, 3, 3, 3, ompgk, &
                3, 3, 1, v3)
!
    do i = 1, 3
        fors(3+i) = v1(i)+v2(i)-v3(i)
    end do
!
    do ne = 1, nno
        coef = -pjacob*en(ne, kp)
        do k = 4, 6
!* ON RETRANCHE LES FORCES D'INERTIE DES FORCES INTERNES:
            fint(k, ne) = fint(k, ne)-coef*fors(k)
        end do
    end do
end subroutine
