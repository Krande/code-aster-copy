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

subroutine cvmcvx(nmat, mater, sig, vin, seuil)
    implicit none
!       VISCOCHABOCHE  :
!          CONVEXE VISCO PLASTIQUE POUR (MATER,SIG,X1,X2,R)  (OPTION 1 )
!                     SEUIL   F    = S   -  AR * R - K
!                                                T           1/2
!                       AVEC  S    = (3/2(D-X1-X2) (D-X1-X2))
!                       ET    D    = SIG - 1/3 TR(SIG) I
!       ----------------------------------------------------------------
!       IN  SIG    :  CONTRAINTE
!       IN  VIN    :  VARIABLES INTERNES = ( X1, X2, P, R, Q, XXI, E3 )
!       IN  NMAT   :  DIMENSION MATER
!       IN  MATER  :  COEFFICIENTS MATERIAU A TEMPERATURE
!       OUT SEUIL  :  SEUIL  ELASTICITE
!       ----------------------------------------------------------------
#include "asterfort/lcdevi.h"
#include "asterfort/lcnrts.h"
    integer(kind=8) :: ndt, ndi, nmat
    real(kind=8) :: sig(6), x1(6), x2(6), dev(6), vin(*)
    real(kind=8) :: difc1, difc2
    real(kind=8) :: ar, r, k, c1d, c2d
    real(kind=8) :: mater(nmat, 2), seuil
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
    common/coed/c1d, c2d
!       ----------------------------------------------------------------
!
! - CALCUL DU PREMIER SEUIL
!
!-----------------------------------------------------------------------
    real(kind=8) :: c1f, c2f
!-----------------------------------------------------------------------
    ar = mater(3, 2)
    k = mater(4, 2)
!        C1D      = MATERD(15,2)
!        C2D      = MATERD(20,2)
    c1f = mater(15, 2)
    c2f = mater(20, 2)
    x1(1:ndt) = vin(1:ndt)
    x2(1:ndt) = vin(ndt+1:ndt+ndt)
!
! --   CAS ANISOTHERME
!
    if (c1d .ne. 0.d0) then
        difc1 = c1f/c1d
        x1(1:ndt) = difc1*x1(1:ndt)
    end if
    if (c2d .ne. 0.d0) then
        difc2 = c2f/c2d
        x2(1:ndt) = difc2*x2(1:ndt)
    end if
!
    r = vin(2*ndt+2)
    call lcdevi(sig, dev)
    dev(1:ndt) = dev(1:ndt)-x1(1:ndt)-x2(1:ndt)
    seuil = lcnrts(dev)-ar*r-k
end subroutine
