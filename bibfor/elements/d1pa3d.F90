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
subroutine d1pa3d(angl, irep, passag)
!.======================================================================
    implicit none
!
!      D1PA3D  -- CALCUL DE LA MATRICE DE PASSAGE DU REPERE
!                 D'ORTHOTROPIE AU REPERE GLOBAL POUR L'INVERSE
!                 DE LA MATRICE DE HOOKE.
!                 CETTE MATRICE EST CONSTRUITE EN ECRIVANT
!                 L'INVARIANCE DE  L'ENERGIE DE DEFORMATION ELASTIQUE
!                 LORS D'UN CHANGEMENT DE REPERE.
!
!   ARGUMENT        E/S  TYPE         ROLE
!    ANLG           IN     R        ANGLES NAUTIQUES DEFINISSANT LE REPERE
!                                   D'ORTHOTROPIE
!    IREP           OUT    I        = 0
!                                     SI LE CHANGEMENT DE REPERE EST
!                                     TRIVIAL (I.E. PASSAG = IDENTITE)
!                                   = 1 SINON
!    PASSAG(6,6)    OUT    R        MATRICE DE PASSAGE DU REPERE
!                                   D'ORTHOTROPIE AU REPERE GLOBAL
!                                   POUR LE TENSEUR D'ELASTICITE
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterfort/matrot.h"
    real(kind=8) :: angl(3), passag(6, 6)
! -----  VARIABLES LOCALES
    real(kind=8) :: p(3, 3)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! ---- INITIALISATIONS
!      ---------------
!-----------------------------------------------------------------------
    integer(kind=8) :: irep
    real(kind=8) :: deux, zero
!-----------------------------------------------------------------------
    zero = 0.0d0
    deux = 2.0d0
    irep = 0
!
    p(:, :) = zero
!
    if (angl(1) .eq. zero .and. angl(2) .eq. zero .and. angl(3) .eq. zero) then
        irep = 0
    else
!       construction de la matrice de passage
        call matrot(angl, p)
        irep = 1
    end if
!
! ---- CONSTRUCTION DE LA MATRICE DE PASSAGE  POUR LE TENSEUR
! ---- D'ELASTICITE (QUI EST DU QUATRIEME ORDRE) DU REPERE
! ---- D'ORTHOTROPIE AU REPERE GLOBAL.
! ---- CETTE MATRICE EST CONSTRUITE EN PARTANT DE LA CONSIDERATION QUE
! ----  (SIGMA_GLOB):(EPSILON_GLOB) = (SIGMA_ORTH):(EPSILON_ORTH)
!       ---------------------------------------------------------
    if (irep .eq. 1) then
!
        passag(1, 1) = p(1, 1)*p(1, 1)
        passag(1, 2) = p(1, 2)*p(1, 2)
        passag(1, 3) = p(1, 3)*p(1, 3)
        passag(1, 4) = deux*p(1, 1)*p(1, 2)
        passag(1, 5) = deux*p(1, 1)*p(1, 3)
        passag(1, 6) = deux*p(1, 2)*p(1, 3)
!
        passag(2, 1) = p(2, 1)*p(2, 1)
        passag(2, 2) = p(2, 2)*p(2, 2)
        passag(2, 3) = p(2, 3)*p(2, 3)
        passag(2, 4) = deux*p(2, 1)*p(2, 2)
        passag(2, 5) = deux*p(2, 1)*p(2, 3)
        passag(2, 6) = deux*p(2, 2)*p(2, 3)
!
        passag(3, 1) = p(3, 1)*p(3, 1)
        passag(3, 2) = p(3, 2)*p(3, 2)
        passag(3, 3) = p(3, 3)*p(3, 3)
        passag(3, 4) = deux*p(3, 1)*p(3, 2)
        passag(3, 5) = deux*p(3, 1)*p(3, 3)
        passag(3, 6) = deux*p(3, 2)*p(3, 3)
!
        passag(4, 1) = p(1, 1)*p(2, 1)
        passag(4, 2) = p(1, 2)*p(2, 2)
        passag(4, 3) = p(1, 3)*p(2, 3)
        passag(4, 4) = (p(1, 1)*p(2, 2)+p(1, 2)*p(2, 1))
        passag(4, 5) = (p(1, 1)*p(2, 3)+p(1, 3)*p(2, 1))
        passag(4, 6) = (p(1, 2)*p(2, 3)+p(1, 3)*p(2, 2))
!
        passag(5, 1) = p(1, 1)*p(3, 1)
        passag(5, 2) = p(1, 2)*p(3, 2)
        passag(5, 3) = p(1, 3)*p(3, 3)
        passag(5, 4) = p(1, 1)*p(3, 2)+p(1, 2)*p(3, 1)
        passag(5, 5) = p(1, 1)*p(3, 3)+p(1, 3)*p(3, 1)
        passag(5, 6) = p(1, 2)*p(3, 3)+p(1, 3)*p(3, 2)
!
        passag(6, 1) = p(2, 1)*p(3, 1)
        passag(6, 2) = p(2, 2)*p(3, 2)
        passag(6, 3) = p(2, 3)*p(3, 3)
        passag(6, 4) = p(2, 1)*p(3, 2)+p(2, 2)*p(3, 1)
        passag(6, 5) = p(2, 1)*p(3, 3)+p(2, 3)*p(3, 1)
        passag(6, 6) = p(2, 2)*p(3, 3)+p(3, 2)*p(2, 3)
!
    end if
end subroutine
