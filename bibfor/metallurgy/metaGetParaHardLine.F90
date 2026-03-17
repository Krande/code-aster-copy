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
subroutine metaGetParaHardLine(poum, fami, kpg, ksp, jvMaterCode, &
                               metaType, nbPhase, &
                               young, coef, h)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/rcvalb.h"
!
    character(len=1), intent(in) :: poum
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp
    integer(kind=8), intent(in) :: jvMaterCode
    integer(kind=8), intent(in) :: metaType, nbPhase
    real(kind=8), intent(in) :: young, coef
    real(kind=8), intent(out) :: h(nbPhase)
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get hardening slope (linear)
!
! --------------------------------------------------------------------------------------------------
!
! In  poum         : '-' or '+' for parameters evaluation (previous or current)
! In  fami         : Gauss family for integration point rule
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  jvMaterCode  : coded material address
! In  metaType     : type of metallurgy
! In  nbPhase      : total number of phasis (cold and hot)
! In  young        : Young modulusure
! In  coef         : coefficient before hardening slope
! Out h            : current hardening slope
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropMaxi = 5
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    character(len=16) :: propName(nbPropMaxi)
    integer(kind=8) :: nbProp, iProp
!
! --------------------------------------------------------------------------------------------------
!
    nbProp = nbPhase
!
    if (metaType .eq. META_STEEL) then
        propName(1) = 'F1_D_SIGM_EPSI'
        propName(2) = 'F2_D_SIGM_EPSI'
        propName(3) = 'F3_D_SIGM_EPSI'
        propName(4) = 'F4_D_SIGM_EPSI'
        propName(5) = 'C_D_SIGM_EPSI'

    elseif (metaType .eq. META_ZIRC) then
        propName(1) = 'F1_D_SIGM_EPSI'
        propName(2) = 'F2_D_SIGM_EPSI'
        propName(3) = 'C_D_SIGM_EPSI'

    else
        ASSERT(ASTER_FALSE)
    end if
!
    call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                ' ', 'META_ECRO_LINE', 0, ' ', [0.d0], &
                nbPhase, propName, propVale, propCode, 2)
    do iProp = 1, nbProp
        h(iProp) = propVale(iProp)
        h(iProp) = coef*h(iProp)*young/(young-h(iProp))
    end do
!
end subroutine
