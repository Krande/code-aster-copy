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
subroutine metaGetParaMixture(poum, fami, kpg, ksp, jvMaterCode, &
                              l_visc, metaType, nbPhase, zalpha, fmix, &
                              sy)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/rcvalb.h"
!
    character(len=1), intent(in) :: poum
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, jvMaterCode
    integer(kind=8), intent(in) :: metaType
    integer(kind=8), intent(in) :: nbPhase
    aster_logical, intent(in) :: l_visc
    real(kind=8), intent(in) :: zalpha
    real(kind=8), intent(out) :: fmix
    real(kind=8), optional, intent(out) :: sy(nbPhase)
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get parameters for mixing law
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
! In  l_visc       : .true. if visco-plasticity
! In  zalpha       : sum of "cold" phasis
! Out fmix         : mixing function
! Out sy           : elasticity yield by phasis
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropMaxi = 6
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    character(len=16) :: propName(nbPropMaxi)
    integer(kind=8) :: nbProp, iProp
!
! --------------------------------------------------------------------------------------------------
!

! - Get mixing function
    nbProp = 1
    propName(1) = 'SY_MELANGE'
    if (l_visc) then
        propName(1) = 'S_VP_MELANGE'
    end if
    fmix = 0.d0
    call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                ' ', 'ELAS_META', 1, 'META', [zalpha], &
                nbProp, propName, propVale, propCode, 0)
    if (propCode(1) .eq. 0) then
        fmix = propVale(1)
    else
        fmix = zalpha
    end if

! - Get elasticity yield by phasis
    if (present(sy)) then
        nbProp = nbPhase
        if (metaType .eq. META_STEEL) then
            propName(1) = 'F1_SY'
            propName(2) = 'F2_SY'
            propName(3) = 'F3_SY'
            propName(4) = 'F4_SY'
            propName(5) = 'C_SY'
            if (l_visc) then
                propName(1) = 'F1_S_VP'
                propName(2) = 'F2_S_VP'
                propName(3) = 'F3_S_VP'
                propName(4) = 'F4_S_VP'
                propName(5) = 'C_S_VP'
            end if
            sy(1:nbProp) = 0.d0

        elseif (metaType .eq. META_ZIRC) then
            propName(1) = 'F1_SY'
            propName(2) = 'F2_SY'
            propName(3) = 'C_SY'
            if (l_visc) then
                propName(1) = 'F1_S_VP'
                propName(2) = 'F2_S_VP'
                propName(3) = 'C_S_VP'
            end if
            sy(1:nbProp) = 0.d0

        else
            ASSERT(ASTER_FALSE)
        end if

        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'ELAS_META', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        do iProp = 1, nbProp
            sy(iProp) = propVale(iProp)
        end do
    end if
!
end subroutine
