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
subroutine metaGetParaAnneal(poum, fami, kpg, ksp, jvMaterCode, &
                             metaType, nbPhaseCold, &
                             annealTheta)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/Metallurgy_type.h"
!
    character(len=1), intent(in) :: poum
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
    integer(kind=8), intent(in) :: metaType, nbPhaseCold
    real(kind=8), intent(out) :: annealTheta(2*nbPhaseCold)
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get parameters for annealing
!
! --------------------------------------------------------------------------------------------------
!
! In  poum         : '-' or '+' for parameters evaluation (previous or current)
! In  fami         : Gauss family for integration point rule
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  jvMaterCode  : coded material address
! In  metaType     : type of metallurgy
! In  nbPhaseCold  : number of cold phases
! Out annealTheta  : parameters for annealing
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropMaxi = 8
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    character(len=16) :: propName(nbPropMaxi)
    integer(kind=8) :: nbProp, iProp
!
! --------------------------------------------------------------------------------------------------
!
    nbProp = 2*(nbPhaseCold-1)

! - Name of parameters
    if (metaType .eq. META_STEEL) then
        propName(1) = 'C_F1_THETA'
        propName(2) = 'C_F2_THETA'
        propName(3) = 'C_F3_THETA'
        propName(4) = 'C_F4_THETA'
        propName(5) = 'F1_C_THETA'
        propName(6) = 'F2_C_THETA'
        propName(7) = 'F3_C_THETA'
        propName(8) = 'F4_C_THETA'
        annealTheta(1:nbProp) = 0.d0
    elseif (metaType .eq. META_ZIRC) then
        propName(1) = 'C_F1_THETA'
        propName(2) = 'C_F2_THETA'
        propName(3) = 'F1_C_THETA'
        propName(4) = 'F2_C_THETA'
        annealTheta(1:nbProp) = 0.d0
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Get parameters
!
    call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                ' ', 'META_RE', 0, ' ', [0.d0], &
                nbProp, propName, propVale, propCode, 2)
    do iProp = 1, nbProp
        annealTheta(iProp) = propVale(iProp)
    end do
!
end subroutine
