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
subroutine sigtopk1(ndim, Cauchy, F, PK1)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/matinv.h"
!
    integer(kind=8), intent(in)                         :: ndim
    real(kind=8), dimension(6), intent(in)    :: Cauchy
    real(kind=8), dimension(3, 3), intent(in)    :: F
    real(kind=8), dimension(3, 3), intent(out)   :: PK1
! --------------------------------------------------------------------------------------------------
! Conversion of the second Piola-Kirshoff stress tensor to the first Piola-Kirshoff stress tensor
! PK1 = J*Sig*F^-T
!
! In ndim: dimension of tthe problem ( 2 or 3)
! In Cauchy : the second Piola-Kirshoff stress tensor (S11, S22, S33, rac2*S12, rac2*S13, rac2*S23)
! In F   : gradient of the deformation
! Out PK1: the first Piola-Kirshoff stress tensor
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8), dimension(3, 3) :: Ft, Ftinv
    real(kind=8) :: J
    real(kind=8), dimension(3, 3) :: SigMat
!
    SigMat(1, 1) = Cauchy(1)
    SigMat(2, 2) = Cauchy(2)
    SigMat(3, 3) = Cauchy(3)
!
    SigMat(1, 2) = Cauchy(4)/rac2
    SigMat(2, 1) = SigMat(1, 2)
!
    if (ndim == 2) then
        SigMat(3, 1:2) = 0.d0
        SigMat(1:2, 3) = 0.d0
    elseif (ndim == 3) then
        SigMat(1, 3) = Cauchy(5)/rac2
        SigMat(3, 1) = SigMat(1, 3)
!
        SigMat(2, 3) = Cauchy(6)/rac2
        SigMat(3, 2) = SigMat(2, 3)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    Ft = transpose(F)
    Ftinv = 0.d0
    call matinv('S', 3, Ft, Ftinv, J)
    if (ndim == 2) Ftinv(3, 3) = 1.d0
!
    PK1 = J*matmul(SigMat, Ftinv)
!
end subroutine
