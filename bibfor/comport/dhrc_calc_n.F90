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

subroutine dhrc_calc_n(eps, vint, b, c, neta1, neta2)
!
! person_in_charge: sebastien.fayolle at edf.fr
!
    implicit none
!
#include "asterfort/r8inir.h"
    real(kind=8), intent(in) :: vint(*), eps(6)
    real(kind=8), intent(in) :: b(6, 2, 2), c(2, 2, 2)
    real(kind=8), intent(out) :: neta1(2), neta2(2)
!
! ----------------------------------------------------------------------
!
!      CALCUL DES FORCES THERMODYNAMIQUES ASSOCIEES A LA PLASTICITE
!      APPELE PAR "SEUGLC"
!
! IN:
!       EPS     : TENSEUR DE DEFORMATIONS
!                 (EXX EYY 2EXY KXX KYY 2KXY)
!       VINT    : VECTEUR DES VARIABLES INTERNES
!                 VINT=(D1,D2,EPSP1X,EPSP1Y,EPSP2X,EPSP2Y)
!       B       : TENSEUR ASSOCIE AUX DEFORMATIONS PLASTIQUES
!       C       : TENSEUR DE RAIDEUR D'ECROUISSAGE
!                 LA PREMIERE COMPOSANTE DE C CORRESPOND AUX GLISSEMENTS
!                 LA DEUXIEME COMPOSANTE DE C CORRESPOND AUX GLISSEMENTS
!                 LA TROISIEME COMPOSANTE DE C CORRESPOND A LA DISTINCTION
!
! OUT:
!       NETA1   : CONTRAINTE ASSOCIEE AU GLISSEMENT SUR LA PARTIE 1 DE
!                 LA PLAQUE
!       NETA2   : CONTRAINTE ASSOCIEE AU GLISSEMENT SUR LA PARTIE 2 DE
!                 LA PLAQUE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: i, k
!
!     INITIALISATION
    call r8inir(2, 0.0d0, neta1, 1)
    call r8inir(2, 0.0d0, neta2, 1)
!
!
    do k = 1, 2
!     CALCUL DE NETA1 ET NETA2
        do i = 1, 2
            neta1(k) = neta1(k)-eps(i)*b(i, k, 1)
            neta2(k) = neta2(k)-eps(i)*b(i, k, 2)
            neta1(k) = neta1(k)-vint(i+2)*c(i, k, 1)
            neta2(k) = neta2(k)-vint(i+4)*c(i, k, 2)
        end do
        do i = 3, 6
            neta1(k) = neta1(k)-eps(i)*b(i, k, 1)
            neta2(k) = neta2(k)-eps(i)*b(i, k, 2)
        end do
    end do
!
end subroutine
