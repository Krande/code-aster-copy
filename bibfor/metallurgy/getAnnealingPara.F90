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
subroutine getAnnealingPara(fami, jvMaterCode, kpg, ksp, &
                            T1, temp, epsqMini, &
                            alpha, tauInf, &
                            lHardMixed, prager, pragerTempEcroIni)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvalb.h"
#include "asterfort/getMixedHardeningPara.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: jvMaterCode
    integer(kind=8), intent(in) :: kpg, ksp
    real(kind=8), intent(in) :: epsqMini, T1, temp
    real(kind=8), intent(out) :: alpha, tauInf
    aster_logical, intent(in) :: lHardMixed
    real(kind=8), intent(out) :: prager, pragerTempEcroIni
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 2
    character(len=16), parameter :: propName(nbProp) = (/"COEF_ECRO", "TAU_INF  "/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    character(len=1), parameter :: poum = '+'
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName = "EPSI"
    real(kind=8) :: paraVale
!
! --------------------------------------------------------------------------------------------------
!
    paraVale = epsqMini
    call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                ' ', 'REST_ECRO', &
                nbPara, paraName, [paraVale], &
                nbProp, propName, propVale, propCode, 2)
    alpha = propVale(1)
    tauInf = propVale(2)
    prager = 0.d0
    pragerTempEcroIni = 0.d0
    if (lHardMixed) then
        call getMixedHardeningPara(fami, jvMaterCode, kpg, ksp, &
                                   T1, prager_=pragerTempEcroIni)
        call getMixedHardeningPara(fami, jvMaterCode, kpg, ksp, &
                                   temp, prager_=prager)
    end if
!
end subroutine
