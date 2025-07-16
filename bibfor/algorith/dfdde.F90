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
subroutine dfdde(eps, endo, ndim, lambda, mu, &
                 dfde)
!
!
    implicit none
#include "asterfort/diago3.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: ndim
    real(kind=8) :: eps(6), lambda, mu
    real(kind=8) :: dfde(6), endo
!
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT ENDO_ORTH_BETON
!     CALCUL DE LA DERIVEE DE LA FORCE THERMODYNAMIQUE(ENDO COMPRESSION)
!     PAR RAPPORT A LA DEFORMATION:DFD/DEPS
!
!     FD=(1-ENDO)(LAMBDA/2*(TR(E).H(-TR(E)))**2+MU*TR(E-**2))-ECROD*ENDO
!     IN  NDIM      : DIMENSION 3(3D) OU 2(2D)
!     IN  EPS        : DEFORMATION
!     IN  ENDO     : ENDOMMAGEMENT DE COMPRESSION
!     IN  LAMBDA   : /
!     IN  MU       : / COEFFICIENTS DE LAME
!     OUT DFDE      : DFD/DEPS
! ----------------------------------------------------------------------
!
    real(kind=8) :: treps, rtemp, rac2
    real(kind=8) :: tu(6), vpe(3)
    real(kind=8) :: valeps(3), veceps(3, 3), phid
    integer(kind=8) :: i, j, k, t(3, 3)
!
    t(1, 1) = 1
    t(1, 2) = 4
    t(1, 3) = 5
    t(2, 1) = 4
    t(2, 2) = 2
    t(2, 3) = 6
    t(3, 1) = 5
    t(3, 2) = 6
    t(3, 3) = 3
!
    rac2 = sqrt(2.d0)
!
!
    phid = 2.d0*(1.d0-endo)
!
    call r8inir(6, 0.d0, dfde, 1)
!
    treps = 0.d0
    treps = eps(1)+eps(2)+eps(3)
!
    if (treps .lt. 0.d0) then
        do i = 1, ndim
            dfde(t(i, i)) = dfde(t(i, i))+phid*lambda*treps
        end do
    end if
!
    call diago3(eps, veceps, valeps)
    call r8inir(3, 0.d0, vpe, 1)
!
    do i = 1, ndim
        if (valeps(i) .lt. 0.d0) then
            vpe(i) = valeps(i)
        else
            vpe(i) = 0.d0
        end if
    end do
!
    call r8inir(6, 0.d0, tu, 1)
    do i = 1, ndim
        do j = i, ndim
            do k = 1, ndim
                tu(t(i, j)) = tu(t(i, j))+veceps(i, k)*vpe(k)*veceps(j, k)
            end do
        end do
    end do
!
    do i = 1, ndim
        do j = i, ndim
            if (j .eq. i) then
                rtemp = 1.d0
            else
                rtemp = rac2
            end if
            dfde(t(i, j)) = dfde(t(i, j))+2.d0*mu*phid*tu(t(i, j))*rtemp
        end do
    end do
!
!
end subroutine
