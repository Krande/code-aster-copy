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
subroutine mctge2(deigy, dydx, direig, eigx, eigy, &
                  edge, outofp)
! ***********************************************************************
!
! OBJECT: COMPUTE THE DERIVATIVE OF A GENERAL ISOTROPIC 2D TENSOR
!
! ----------------------------------------------------------------------
!
!     LOI DE COMPORTEMENT DE MOHR-COULOMB
!
! IN  DEIGY   : D Y_PRIN / D X_PRIN (3,3)-MATRIX
! IN  DIREIG  : DIRECTIONS PRINCIPALES DE Y_PRIN:
!                  DIREIG_1 = DIREIG(I=1-3,1)
!                  DIREIG_2 = DIREIG(I=1-3,2)
!                  DIREIG_3 = DIREIG(I=1-3,3)
! IN  EIGX    : VALEURS PRINCIPALES DE X (3)
! IN  EIGY    : VALEURS PRINCIPALES DE Y (3)
! IN  EDGE    : Y-A-T-IL DEUX  MECANISMES ACTIFS
! IN  OUTOFP  : COMPOSANTE HORS PLAN = | TRUE :  D_PLAN
!                                      | FALSE:  C_PLAN
!
! OUT DYDX    : MATRICE TANGENTE COHERENTE REACTUALISEE
!
! ***********************************************************************
    implicit none
! ======================================================================
!
!
#include "asterf_types.h"
!
    real(kind=8) :: deigy(3, 3)
    real(kind=8) :: dydx(6, 6)
    real(kind=8) :: direig(3, 3)
    real(kind=8) :: eigx(3)
    real(kind=8) :: eigy(3)
    real(kind=8) :: edge
    aster_logical :: outofp
!
! Declaration of integer type variables
    integer(kind=8) :: i, j, mcomp, mdim, ndim
!
    parameter(mcomp=4, mdim=3, ndim=2)
!
    real(kind=8) :: r1, r2, r3, r4, rp5, eigpr3(mcomp), foid(mcomp, mcomp)
    real(kind=8) :: sopid(mcomp), a1, sqr, r0, eigprj(mcomp, ndim)
    data r0, r1, r2, r3, r4, rp5, sqr/&
     &    0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 0.5d0,&
     &    1.4142135623730951d0/
!
    dydx(:, :) = r0
    foid(:, :) = r0
    eigprj(:, :) = r0
    sopid(:) = r0
    eigpr3(:) = r0
!
    do i = 1, ndim
        foid(i, i) = r1
        sopid(i) = r1
    end do
    foid(mcomp, mcomp) = rp5
    eigpr3(mdim) = r1
!
! Calculation of the eigenvectors EIGPRJ_1(6) EIGPRJ_2(6) EIGPRJ_3(6)
! from the eigendirections DIREIG_1(3) DIREIG_2(3) DIREIG_3(3)
! EIGPRJ_A_ij = DIREIG_A_i x DIREIG_A_j etc.
    do i = 1, mdim
        eigprj(i, 1) = direig(i, 1)*direig(i, 1)
        eigprj(i, 2) = direig(i, 2)*direig(i, 2)
    end do
    eigprj(mcomp, 1) = direig(1, 1)*direig(2, 1)
    eigprj(mcomp, 2) = direig(1, 2)*direig(2, 2)
!
    if (edge .eq. r1) then
!
!
! Derivative dY/dX for repeated in-plane eigenvalues of X
! -------------------------------------------------------
! In-plane component
        do i = 1, mcomp
            do j = 1, mcomp
                dydx(i, j) = (deigy(1, 1)-deigy(1, 2))*foid(i, j)+ &
                             deigy(1, 2)*sopid(i)*sopid(j)
            end do
        end do
!
! out-of-plane components required
        if (outofp) then
            do i = 1, mcomp
                do j = 1, mcomp
                    if (i .eq. mdim .or. j .eq. mdim) then
                        dydx(i, j) = deigy(1, 3)*sopid(i)*eigpr3(j)+ &
                                     deigy(3, 1)*eigpr3(i)*sopid(j)+ &
                                     deigy(3, 3)*eigpr3(i)*eigpr3(j)
                    end if
                end do
            end do
        end if
    else
!
! Derivative dY/dX for distinct in-plane eigenvalues of X
! -------------------------------------------------------
! Assemble in-plane DYDX
        a1 = (eigy(1)-eigy(2))/(eigx(1)-eigx(2))
        do i = 1, mcomp
            do j = 1, mcomp
                dydx(i, j) = a1*(foid(i, j)- &
                                 eigprj(i, 1)*eigprj(j, 1)- &
                                 eigprj(i, 2)*eigprj(j, 2))+ &
                             deigy(1, 1)*eigprj(i, 1)*eigprj(j, 1)+ &
                             deigy(1, 2)*eigprj(i, 1)*eigprj(j, 2)+ &
                             deigy(2, 1)*eigprj(i, 2)*eigprj(j, 1)+ &
                             deigy(2, 2)*eigprj(i, 2)*eigprj(j, 2)
            end do
        end do
!
! out-of-plane components required
        if (outofp) then
            do i = 1, mcomp
                do j = 1, mcomp
                    if (i .eq. mdim .or. j .eq. mdim) then
                        dydx(i, j) = deigy(1, 3)*eigprj(i, 1)*eigpr3(j)+ &
                                     deigy(2, 3)*eigprj(i, 2)*eigpr3(j)+ &
                                     deigy(3, 1)*eigpr3(i)*eigprj(j, 1)+ &
                                     deigy(3, 2)*eigpr3(i)*eigprj(j, 2)+ &
                                     deigy(3, 3)*eigpr3(i)*eigpr3(j)
                    end if
                end do
            end do
        end if
    end if
!
! Projection to the so-called "tensor base"
    do i = 1, mcomp
        dydx(mcomp, i) = sqr*dydx(mcomp, i)
        dydx(i, mcomp) = sqr*dydx(i, mcomp)
    end do
!
end subroutine
