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
subroutine dfddd(eps, endo, ndim, lambda, mu, &
                 ecrod, dfd)
!
!
    implicit none
#include "asterfort/diago3.h"
#include "asterfort/r8inir.h"
    integer(kind=8) :: ndim
    real(kind=8) :: eps(6), lambda, mu
    real(kind=8) :: dfd, endo, ecrod
!
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT ENDO_ORTH_BETON
!     CALCUL DE LA DERIVEE DE LA FORCE THERMODYNAMIQUE(ENDO COMPRESSION)
!     PAR RAPPORT A L ENDOMMAGEMENT DE COMPRESSION:DFD/DD
!
!     FD=(1-ENDO)(LAMBDA/2*(TR(E).H(-TR(E)))**2+MU*TR(E-**2))-ECROD*ENDO
!     IN  NDIM      : DIMENSION 3(3D) OU 2(2D)
!     IN  EPS      : DEFORMATION
!     IN  ENDO     : ENDOMMAGEMENT DE COMPRESSION
!     IN  LAMBDA   : /
!     IN  MU       : / COEFFICIENTS DE LAME
!     IN  ECROD    : PARAMETRE D ECROUISSAGE ENDO COMPRESSION
!     OUT DFD      : DFD/DD
! ----------------------------------------------------------------------
!
!
    real(kind=8) :: treps
    real(kind=8) :: vpe(3)
    real(kind=8) :: valeps(3), veceps(3, 3), phid
    integer(kind=8) :: i, t(3, 3)
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
    phid = -2.d0
!
    dfd = 0.d0
    treps = 0.d0
    do i = 1, ndim
        treps = treps+eps(t(i, i))
    end do
    if (treps .lt. 0.d0) then
        dfd = dfd+phid*lambda/2.d0*treps**2
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
    dfd = dfd+mu*phid*(vpe(1)**2+vpe(2)**2+vpe(3)**2)-2.d0*ecrod
!
end subroutine
