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

subroutine irrfss(sig, ddfdds)
    implicit none
#include "asterfort/lcdevi.h"
#include "asterfort/lcnrts.h"
#include "asterfort/lcprte.h"
    real(kind=8) :: sig(6), ddfdds(6, 6)
!       3D / DP / CP
!       DERIVEE / S / S DE LA FONCTION SEUIL A (SIG ) DONNES
!
!       IN  SIG    :  TENSEUR CONTRAINTE
!                                                     T
!       OUT D2FD2S :  DDFDDS = 1/S ( 3/2 ID - DFDS DFDS )
!                     DFDS   = 3 D / 2 S
!                                            T
!                     ID     = I4 - 1/3 I2 I2
!                                           T           1/2
!                     S      = (3/2(D-X1-X2) (D-X1-X2))
!                     D      = SIG - 1/3 TR(SIG) I
!       ----------------------------------------------------------------
    real(kind=8) :: id(6, 6), d23, d13, zero, un, s
    parameter(d23=.66666666666666d0)
    parameter(d13=-.33333333333333d0)
    parameter(zero=0.d0)
    parameter(un=1.d0)
    real(kind=8) :: dev(6), dfds(6), dfds2(6, 6)
    integer(kind=8) :: ndt, ndi
    common/tdim/ndt, ndi
    data id/d23, d13, d13, zero, zero, zero,&
     &                  d13, d23, d13, zero, zero, zero,&
     &                  d13, d13, d23, zero, zero, zero,&
     &                  zero, zero, zero, un, zero, zero,&
     &                  zero, zero, zero, zero, un, zero,&
     &                  zero, zero, zero, zero, zero, un/
!
!
!
    call lcdevi(sig, dev)
    s = lcnrts(dev)
    if (s .eq. 0.d0) then
        ddfdds(:, :) = 0.d0
    else
        dfds(1:ndt) = (1.5d0/s)*dev(1:ndt)
        call lcprte(dfds, dfds, dfds2)
        ddfdds(1:ndt, 1:ndt) = 1.5d0*id(1:ndt, 1:ndt)
        ddfdds(1:ndt, 1:ndt) = ddfdds(1:ndt, 1:ndt)-dfds2(1:ndt, 1:ndt)
        ddfdds(1:ndt, 1:ndt) = (1.d0/s)*ddfdds(1:ndt, 1:ndt)
    end if
!
end subroutine
