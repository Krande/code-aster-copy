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
subroutine getAnnealingParaOnCell(fami, jvMaterCode, &
                                  T1, T2, &
                                  epsqMini, xcinMini)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvalb.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: jvMaterCode
    real(kind=8), intent(out) :: T1, T2, epsqMini, xcinMini
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: kpg = 1, kspg = 1
    character(len=1), parameter :: poum = '+'
    integer(kind=8), parameter :: nbProp = 4
    integer(kind=8) :: propCode(nbProp)
    character(len=16), parameter :: propName(nbProp) = &
                                    (/'TEMP_MINI', 'TEMP_MAXI', 'EPSQ_MINI', 'XCIN_MINI'/)
    real(kind=8) :: propVale(nbProp)
!
! --------------------------------------------------------------------------------------------------
!
    call rcvalb(fami, kpg, kspg, poum, &
                jvMaterCode, ' ', 'REST_ECRO', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 2)
    T1 = propVale(1)
    T2 = propVale(2)
    epsqMini = propVale(3)
    xcinMini = propVale(4)
!
end subroutine
