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
subroutine lkijpl(nmat, mater, sigf, nr, drdy, &
                  dsde)
    implicit none
! person_in_charge: alexandre.foucault at edf.fr
!       ----------------------------------------------------------------
!       MATRICE COHERENTE DE LETK A T+DT
!       IN  NMAT   :  DIMENSION MATER
!           MATER  :  COEFFICIENTS MATERIAU
!           SIGF   :  ETAT DE CONTRAINTES A T+DT
!           NR     :  DIMENSION MATRICE JACOBIENNE
!           DRDY   :  MATRICE JACOBIENNE (NR*NR)
!       OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT = DSIG/DEPS
!       ----------------------------------------------------------------
#include "asterc/r8prem.h"
#include "asterfort/lcopli.h"
#include "asterfort/mgauss.h"
#include "asterfort/prmama.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: nmat, nr
    real(kind=8) :: dsde(6, 6), mater(nmat, 2)
    real(kind=8) :: drdy(nr, nr), sigf(6)
!
    integer(kind=8) :: i, j, iret, ier, ndt, ndi
    real(kind=8) :: jss(6, 6), jsz(6, 3), jzs(3, 6), jzz(3, 3)
    real(kind=8) :: hook(6, 6), hooknl(6, 6), i1, coefnl
    real(kind=8) :: patm, nelas, invjzz(3, 3), j3x6(3, 6)
    real(kind=8) :: det, j6x6(6, 6), dijaco(6, 6), invdij(6, 6)
    real(kind=8) :: maxi, mini, mue, mu
!
!       --------------------------------------------------------------
    common/tdim/ndt, ndi
!       --------------------------------------------------------------
! === =================================================================
! --- INITIALISATION MATRICES A ZERO
! === =================================================================
    jss(:, :) = 0.d0
    call r8inir(18, 0.d0, jsz, 1)
    call r8inir(18, 0.d0, jzs, 1)
    call r8inir(9, 0.d0, jzz, 1)
! === =================================================================
! --- RECHERCHE DU MAXIMUM DE DRDY
! === =================================================================
    maxi = 0.d0
    do i = 1, nr
        do j = 1, nr
            if (abs(drdy(i, j)) .gt. maxi) maxi = abs(drdy(i, j))
        end do
    end do
! === =================================================================
! --- DIMENSIONNEMENT A R8PREM
! === =================================================================
    mini = r8prem()*maxi
    do i = 1, nr
        do j = 1, nr
            if (abs(drdy(i, j)) .lt. mini) drdy(i, j) = 0.d0
        end do
    end do
!
! === =================================================================
! --- SEPARATION DES TERMES DU JACOBIEN
! === =================================================================
    do i = 1, ndt
        do j = 1, ndt
            jss(i, j) = drdy(i, j)
        end do
    end do
!
    do i = 1, 3
        do j = 1, ndt
            jsz(j, i) = drdy(j, ndt+i)
            jzs(i, j) = drdy(ndt+i, j)
        end do
    end do
!
    do i = 1, 3
        do j = 1, 3
            jzz(i, j) = drdy(ndt+i, ndt+j)
        end do
    end do
! === =================================================================
! --- CONSTRUCTION TENSEUR RIGIDITE ELASTIQUE A T+DT
! === =================================================================
    call lcopli('ISOTROPE', '3D      ', mater, hook)
! --- PRISE EN COMPTE DU TERME NON LINEAIRE E(I1+) = E0*(I1+/PA)**NE
! --- AVEC I1 = -TRACE(SIGMA), CAR EQUATIONS DU MODELE LETK
! --- SONT EXPRIMEES EN CONVENTION MECANIQUE DES SOLS
    i1 = -(sigf(1)+sigf(2)+sigf(3))
    patm = mater(1, 2)
    nelas = mater(2, 2)
    coefnl = (i1/(3.d0*patm))**nelas
!
    mue = mater(4, 1)
    mu = -mue*coefnl
!
! === =================================================================
! --- MISE A L'ECHELLE DU NUMERATEUR DR(1:6)/DEPS
! === =================================================================
    coefnl = coefnl/mu
!
    hooknl(1:ndt, 1:ndt) = coefnl*hook(1:ndt, 1:ndt)
! === =================================================================
! --- CONSTRUCTION TENSEUR CONSTITUTIF TANGENT DSDE
! === =================================================================
! --- INVERSION DU TERME JZZ
    call r8inir(9, 0.d0, invjzz, 1)
    do i = 1, 3
        invjzz(i, i) = 1.d0
    end do
!
    call mgauss('NCVP', jzz, invjzz, 3, 3, &
                3, det, iret)
    if (iret .gt. 0) call r8inir(9, 0.d0, invjzz, 1)
!
! --- PRODUIT DU TERME (JZZ)^-1*JZS = J3X6
    call prmama(1, invjzz, 3, 3, 3, &
                jzs, 3, 3, ndt, j3x6, &
                3, 3, ndt, ier)
    if (ier .gt. 0) write (6, *) 'ECHEC AVEC PRMAMA 1'
!
! --- PRODUIT DU TERME JSZ*(JZZ)^-1*JZS = JSZ*J3*6 = J6X6
    call prmama(1, jsz, 6, ndt, 3, &
                j3x6, 3, 3, ndt, j6x6, &
                6, ndt, ndt, ier)
    if (ier .gt. 0) write (6, *) 'ECHEC AVEC PRMAMA 2'
!
! --- DIFFERENCE DE MATRICE (JSS - J6X6) = DIJACO
    dijaco(1:ndt, 1:ndt) = jss(1:ndt, 1:ndt)-j6x6(1:ndt, 1:ndt)
!
! --- INVERSION DU TERME (DIJACO)^-1 = INVDIJ
    invdij(:, :) = 0.d0
    do i = 1, ndt
        invdij(i, i) = 1.d0
    end do
    call mgauss('NCVP', dijaco, invdij, 6, ndt, &
                ndt, det, iret)
    if (iret .gt. 1) then
        dsde(1:ndt, 1:ndt) = hook(1:ndt, 1:ndt)
    end if
!
! --- CONSTRUCTION DSDE = INVDIJ*HOOKNL
    dsde(1:ndt, 1:ndt) = matmul(invdij(1:ndt, 1:ndt), hooknl(1:ndt, 1:ndt))
!
end subroutine
