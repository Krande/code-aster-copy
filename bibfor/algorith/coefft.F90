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
subroutine coefft(cothe, coeff, dcothe, dcoeff, x, &
                  dtime, coeft, nmat, coel)
    implicit none
!       INTEGRATION DE LOIS DE COMPORTEMENT ELASTO-VISCOPLASTIQUE PAR
!       UNE METHODE DE RUNGE KUTTA AVEC REDECOUPAGE AUTOMATIQUE DU PAS
!       DE TEMPS : CALCUL DES PARAMETRES MATERIAU A UN INSTANT DONNE
!       ---------------------------------------------------------------
!       IN  COTHE  :  COEFFICIENTS MATERIAU ELASTIQUE A T
!           COEFF  :  COEFFICIENTS MATERIAU INELASTIQUE A T
!           DCOTHE :  INTERVALLE COEFFICIENTS MATERIAU ELAST POUR DT
!           DCOEFF :  INTERVALLE COEFFICIENTS MATERIAU INELAST POUR DT
!           X      :  INSTANT COURANT
!           DTIME  :  INTERVALLE DE TEMPS
!           NMAT   :  NOMNRE MAXI DE COEF MATERIAU
!       OUT COEFT  :  COEFFICIENTS MATERIAU INELASTIQUE A T+DT
!C          COEL   :  COEFFICIENTS  ELASTIQUES ELASTIQUE A T+DT
!       ---------------------------------------------------------------
#include "asterfort/r8inir.h"
    integer(kind=8) :: nmat, i, ncoe, ncoel
    real(kind=8) :: cothe(nmat), dcothe(nmat), coel(nmat)
    real(kind=8) :: hsdt, dtime, x
    real(kind=8) :: coeff(nmat), dcoeff(nmat), coeft(nmat)
!
    hsdt = x/dtime
!
    call r8inir(nmat, 0.d0, coel, 1)
!
    if (cothe(nmat) .eq. 0) then
        do i = 1, 3
            coel(i) = cothe(i)+hsdt*dcothe(i)
        end do
        coel(nmat) = 0.d0
    else if (cothe(nmat) .eq. 1) then
!          MATERIAU ISOTROPE, OU ANISOTROPE, MATRICE DE HOOKE
        coel(nmat) = 1.d0
        ncoel = 75
        do i = 1, ncoel
            coel(i) = cothe(i)+hsdt*dcothe(i)
        end do
    end if
!
!       POUR GAGNER DU TEMPS CPU
    if (coeff(nmat) .eq. 0) then
        ncoe = nmat
    else
        ncoe = nint(coeff(nmat))
    end if
!
    do i = 1, ncoe
        coeft(i) = coeff(i)+hsdt*dcoeff(i)
    end do
!
end subroutine
