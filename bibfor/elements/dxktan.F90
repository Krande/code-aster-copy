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
subroutine dxktan(delas, mp1, mp2, nbackn, ncrit, &
                  dcc1, dcc2, dsidep)
    implicit none
!     REALISE LE CALCUL DE LA MATRICE TANGENTE DANS LE CAS DE LA LOI
!     DE COMPORTEMENT GLRC
!
! IN  DELAS  : MATRICE ELASTIQUE EN MEMBRANE, FLEXION ET COUPLAGE
! IN  MP1    : MOMENTS LIMITES ELASTIQUES EN FLEXION POSITIVE
! IN  MP2    : MOMENTS LIMITES ELASTIQUES EN FLEXION NEGATIVE
! IN  NBACKN : MOMENT DE RAPPEL
! IN  NCRIT  : TYPE DU CRITERE DE PLASTICITE
! IN  DCC1   : MATRICE ELASTIQUE + CONSTANTES DE PRAGER (FLEXION +)
! IN  DCC2   : MATRICE ELASTIQUE + CONSTANTES DE PRAGER (FLEXION -)
!
! OUT DSIDEP : MATRICE TANGENTE
! ----------------------------------------------------------------------
!
#include "asterfort/dfplas.h"
#include "asterfort/dxprd1.h"
#include "asterfort/dxprd2.h"
#include "asterfort/lcprte.h"
#include "asterfort/pmavec.h"
#include "blas/dcopy.h"
    integer(kind=8) :: ncrit, i, j, n, nd
!
    real(kind=8) :: mp1(3), mp2(3), nbackn(6)
    real(kind=8) :: dfpla1(6)
    real(kind=8) :: dfpla2(6)
    real(kind=8) :: delas(6, 6), dsidep(6, 6)
    real(kind=8) :: mat(6, 6), mata(6, 6), matb(6, 6), matc(6, 6), matd(6, 6)
    real(kind=8) :: dcc1(3, 3), dcc2(3, 3), dc1(6, 6), dc2(6, 6)
    real(kind=8) :: vect(6)
    real(kind=8) :: scal, scala, scalb
    blas_int :: b_incx, b_incy, b_n
    common/tdim/n, nd
!
!     INITIALISATION
    n = 6
!
    do i = 1, 6
        dfpla1(i) = 0.d0
        dfpla2(i) = 0.d0
        do j = 1, 6
            dc1(i, j) = 0.d0
            dc2(i, j) = 0.d0
        end do
    end do
!
    b_n = to_blas_int(36)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, delas, b_incx, dc1, b_incy)
    b_n = to_blas_int(36)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, delas, b_incx, dc2, b_incy)
!
    do i = 1, 3
        do j = 1, 3
            dc1(i+3, j+3) = dcc1(i, j)
            dc2(i+3, j+3) = dcc2(i, j)
        end do
    end do
!
    if (ncrit .eq. 0) then
!     CAS ELASTIQUE
        b_n = to_blas_int(36)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, delas, b_incx, dsidep, b_incy)
!
    else if (ncrit .eq. 1) then
        call dfplas(nbackn(4), mp1, dfpla1(4))
        call lcprte(dfpla1, dfpla1, mata)
        matb = matmul(delas, mata)
        matc = matmul(matb, delas)
        call pmavec('ZERO', 6, dc1, dfpla1, vect)
        scal = dot_product(dfpla1(1:n), vect(1:n))
!
        do i = 1, 6
            do j = 1, 6
                dsidep(i, j) = delas(i, j)-matc(i, j)/scal
            end do
        end do
!
    else if (ncrit .eq. 2) then
        call dfplas(nbackn(4), mp2, dfpla2(4))
        call lcprte(dfpla2, dfpla2, mata)
        matb = matmul(delas, mata)
        matc = matmul(matb, delas)
        call pmavec('ZERO', 6, dc2, dfpla2, vect)
        scal = dot_product(dfpla2(1:n), vect(1:n))
!
        do i = 1, 6
            do j = 1, 6
                dsidep(i, j) = delas(i, j)-matc(i, j)/scal
            end do
        end do
    else if (ncrit .eq. 12) then
!
!     NUMERATEUR
        call dfplas(nbackn(4), mp1, dfpla1(4))
        call dfplas(nbackn(4), mp2, dfpla2(4))
        call dxprd1(dfpla1, dfpla2, dc2, dfpla2, dfpla1, &
                    mata)
        call dxprd1(dfpla1, dfpla1, dc2, dfpla2, dfpla2, &
                    matb)
        call dxprd1(dfpla2, dfpla1, dc1, dfpla1, dfpla2, &
                    matc)
        call dxprd1(dfpla2, dfpla2, dc1, dfpla1, dfpla1, &
                    matd)
!     DENOMINATEUR
        call dxprd2(dfpla1, dc1, dfpla1, dfpla2, dc2, &
                    dfpla2, scala)
        call dxprd2(dfpla2, dc1, dfpla1, dfpla1, dc2, &
                    dfpla2, scalb)
!
        scal = scala*scalb
!
        do i = 1, 6
            do j = 1, 6
                mat(i, j) = (mata(i, j)-matb(i, j)+matc(i, j)-matd(i, j))/scal
            end do
        end do
!
        mata = matmul(delas, mat)
        matb = matmul(mata, delas)
!
        do i = 1, 6
            do j = 1, 6
                dsidep(i, j) = delas(i, j)-matb(i, j)
            end do
        end do
!
    end if
end subroutine
