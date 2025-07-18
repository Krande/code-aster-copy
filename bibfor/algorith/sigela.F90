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
subroutine sigela(typmod, ndim, e, nu, epse, &
                  sigel)
!
    implicit none
#include "asterfort/bptobg.h"
#include "asterfort/jacobi.h"
#include "asterfort/r8inir.h"
    character(len=8) :: typmod(1)
    integer(kind=8) :: ndim
    real(kind=8) :: epse(6), e, nu
    real(kind=8) :: sigel(6)
! ----------------------------------------------------------------------
!  CALCUL DES CONTRAINTES ELASTIQUES A PARTIR DES DEFORMATIONS
!   ELASTIQUES POUR UN MATERIAU ISOTROPE
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  EPSE    : DEFORMATIONS ELASTIQUES
! IN  E       : MODULE ELASTIQUE DE YOUNG
! IN  NU      : COEFFICIENT DE POISSON
! IN TYPMOD   : TYPE DE MODELISATION (C_PLAN...)
!
! OUT SIGEL   : CONTRAINTES ELASTIQUES
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ndimsi, nperm, nitjac, trij, ordrej, k
    real(kind=8) :: epsep(3), vecpe(3, 3), sigelp(3)
    real(kind=8) :: tol, toldyn, tr(6), tu(6), jacaux(3)
    real(kind=8) :: rac2, coplan, lambda, deuxmu
!
! ======================================================================
!
!
!--------------------------------------------------------
!                            INITIALISATION
!--------------------------------------------------------
!
    ndimsi = 2*ndim
    rac2 = sqrt(2.d0)
!
    lambda = e*nu/(1.d0+nu)/(1.d0-2.d0*nu)
    deuxmu = e/(1.d0+nu)
!
    if (typmod(1) .eq. 'C_PLAN  ') then
        coplan = -nu/(1.d0-nu)
        epse(3) = coplan*(epse(1)+epse(2))
    end if
!
    do k = 4, ndimsi
        epse(k) = epse(k)/rac2
    end do
!
!--------------------------------------------------------
!  -   ON PASSE DANS LE REPERE PROPRE DE EPS
!--------------------------------------------------------
!
    nperm = 12
    tol = 1.d-10
    toldyn = 1.d-2
!       MATRICE  TR = (XX XY XZ YY YZ ZZ) POUR JACOBI)
    tr(1) = epse(1)
    tr(2) = epse(4)
    tr(3) = epse(5)
    tr(4) = epse(2)
    tr(5) = epse(6)
    tr(6) = epse(3)
!     MATRICE UNITE = (1 0 0 1 0 1) (POUR JACOBI)
    tu(1) = 1.d0
    tu(2) = 0.d0
    tu(3) = 0.d0
    tu(4) = 1.d0
    tu(5) = 0.d0
    tu(6) = 1.d0
    trij = 2
    ordrej = 2
!
    call jacobi(3, nperm, tol, toldyn, tr, &
                tu, vecpe, epsep, jacaux, nitjac, &
                trij, ordrej)
!
!
!----------------------------------------------------------------
!     CALCUL DES CONTRAINTES ELASTIQUES (REPERE PRINCIPAL)
!----------------------------------------------------------------
    do k = 1, 3
        sigelp(k) = lambda*(epsep(1)+epsep(2)+epsep(3))
    end do
    do k = 1, 3
        sigelp(k) = sigelp(k)+deuxmu*epsep(k)
    end do
!
!------------------------------------------------------------------
!     ON PASSE DANS LE REPERE INITIAL LES CONTRAINTES ELASTIQUES
!------------------------------------------------------------------
    call r8inir(6, 0.d0, sigel, 1)
    tr(1) = sigelp(1)
    tr(2) = sigelp(2)
    tr(3) = sigelp(3)
    tr(4) = 0.d0
    tr(5) = 0.d0
    tr(6) = 0.d0
    call bptobg(tr, sigel, vecpe)
    do k = 4, ndimsi
        sigel(k) = rac2*sigel(k)
    end do
end subroutine
