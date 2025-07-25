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
subroutine lccong(nr, itmax, toler, iter, r, &
                  rini, yd, dy, irtet)
!     CONTROLE DE LA CONVERGENCE DU NEWTON LOCAL DE LETK
!                   - CONTROLE DU NOMBRE D ITERATIONS
!                   - CONTROLE DE LA PRECISION DE CONVERGENCEC
!     ----------------------------------------------------------------
!     IN  ITMAX  :  NB MAXI D ITERATIONS LOCALES
!         TOLER  :  TOLERANCE A CONVERGENCE
!         ITER   :  NUMERO ITERATION COURANTE
!         NR     :  DIMENSION R
!         R      :  RESIDU DU SYSTEME NL A L'ITERATION COURANTE
!         RINI   :  RESIDU DU SYSTEME NL A LA 1ERE ITERATION
!         YD     :  SOLUTION A DEBUT DU PAS DE TEMPS
!         DY     :  INCREMENT DE SOLUTION
!
!     OUT IRET = 0:  CONVERGENCE
!         IRET = 1:  ITERATION SUIVANTE
!         IRET = 2:  RE-INTEGRATION
!         IRET = 3:  REDECOUPAGE DU PAS DE TEMPS
!     ----------------------------------------------------------------
    implicit none
#include "asterc/r8prem.h"
    integer(kind=8) :: nr, itmax, iter, irtet, i
    real(kind=8) :: toler, r(nr), e1, e2, e1ini, e2ini, errr(2), rini(*)
    real(kind=8) :: yd(*), dy(*), err
!     ----------------------------------------------------------------
! === ==================================================================
! --- CALCUL DE LA NORME DE RINI ET DE R(Y)
! === ==================================================================
!
    e1 = 0.d0
    e1ini = 0.d0
    do i = 1, 6
        e1 = max(e1, abs(r(i)))
        e1ini = max(e1ini, abs(rini(i)))
    end do
!     R8PREM CAR R HOMOGENE A DES DEFORMATIONS
    errr(1) = e1
    if (e1ini .gt. r8prem()) then
        errr(1) = e1/e1ini
    end if
!
    e2 = 0.d0
    e2ini = 0.d0
    do i = 7, nr
        e2 = max(e2, abs(r(i)))
        e2ini = max(e2ini, abs(yd(i)+dy(i)))
    end do
!
    errr(2) = e2
    if (e2ini .gt. r8prem()) then
        errr(2) = e2/e2ini
    end if
!
!     MAX DES 6 PREMIERS TERMES ET DES SUIVANTS
    err = max(errr(1), errr(2))
!
! === =================================================================
! --- TEST DE CONVERGENCE PAR RAPPORT A TOLER
! === =================================================================
    if (err .lt. toler) then
        irtet = 0
        goto 999
    end if
!
! === ==================================================================
! --- SI NON CONVERGENCE: TEST DU N°ITERATION
! === ==================================================================
    if (iter .lt. itmax) then
        irtet = 1
    else
        irtet = 3
    end if
!
999 continue
!
end subroutine
