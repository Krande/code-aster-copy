! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine getMixedHardeningPara(fami, jvMaterCode, kpg, ksp, &
                                 temp, young_, prager_, dsde_, sigy_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvalb.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: jvMaterCode
    integer(kind=8), intent(in) :: kpg, ksp
    real(kind=8), intent(in) :: temp
    real(kind=8), optional, intent(out) :: young_, prager_, dsde_, sigy_
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbProp
    character(len=16) :: propName(3)
    real(kind=8) :: propVale(3)
    integer(kind=8) :: propCode(3)
    character(len=1), parameter :: poum = '+'
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName = "TEMP"
    real(kind=8) :: paraVale
    real(kind=8) :: young, prager, dsde, sigy
!
! --------------------------------------------------------------------------------------------------
!
    paraVale = temp

! - Get elasticity
    nbProp = 1
    propName(1) = "E"
    call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                ' ', 'ELAS', nbPara, paraName, [paraVale], &
                nbProp, propName, propVale, propCode, 2)
    young = propVale(1)

! - Get Prager
    nbProp = 1
    propName(1) = "C"
    call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                ' ', 'PRAGER', nbPara, paraName, [paraVale], &
                nbProp, propName, propVale, propCode, 2)
    prager = propVale(1)

! - Get kinematic parameters
    nbProp = 2
    propName(1) = 'D_SIGM_EPSI'
    propName(2) = 'SY'
    call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                ' ', 'ECRO_LINE', nbPara, paraName, [paraVale], &
                nbProp, propName, propVale, propCode, 2)
    dsde = propVale(1)
    sigy = propVale(2)
!
    if (present(young_)) then
        young_ = young
    end if
    if (present(dsde_)) then
        dsde_ = dsde
    end if
    if (present(sigy_)) then
        sigy_ = sigy
    end if
    if (present(prager_)) then
        prager_ = prager
    end if
!
end subroutine
