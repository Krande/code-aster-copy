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
subroutine metaGetParaVisc(poum, fami, kpg, ksp, jvMaterCode, &
                           metaType, nbPhase, &
                           eta, n, unsurn, &
                           c, m)
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
    integer(kind=8), intent(in) :: metaType
    integer(kind=8), intent(in) :: nbPhase
    real(kind=8), optional, intent(out) :: eta(nbPhase), n(nbPhase), unsurn(nbPhase)
    real(kind=8), optional, intent(out) :: c(nbPhase), m(nbPhase)
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get parameters for viscosity
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
! Out eta          : viscosity parameter - eta
! Out n            : viscosity parameter - n
! Out unsurn       : viscosity parameter - 1/n
! Out c            : viscosity parameter - C
! Out m            : viscosity parameter - m
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropMaxi = 5
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    character(len=8) :: propName(nbPropMaxi)
    integer(kind=8) :: nbProp, iProp
!
! --------------------------------------------------------------------------------------------------
!
    nbProp = nbPhase

! - Name of parameters
    if (metaType .eq. META_STEEL) then
        if (present(eta)) then
            propName(1) = 'F1_ETA'
            propName(2) = 'F2_ETA'
            propName(3) = 'F3_ETA'
            propName(4) = 'F4_ETA'
            propName(5) = 'C_ETA'
            eta(1:nbProp) = 0.d0
        end if
    elseif (metaType .eq. META_ZIRC) then
        if (present(eta)) then
            propName(1) = 'F1_ETA'
            propName(2) = 'F2_ETA'
            propName(3) = 'C_ETA'
            eta(1:nbProp) = 0.d0
        end if
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get parameters
    if (present(eta)) then
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'META_VISC', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        do iProp = 1, nbProp
            eta(iProp) = propVale(iProp)
        end do
    end if

! - Name of parameters
    if (metaType .eq. META_STEEL) then
        if (present(n)) then
            propName(1) = 'F1_N'
            propName(2) = 'F2_N'
            propName(3) = 'F3_N'
            propName(4) = 'F4_N'
            propName(5) = 'C_N'
            n(1:nbProp) = 20.d0
            unsurn(1:nbProp) = 1.d0
        end if
    elseif (metaType .eq. META_ZIRC) then
        if (present(n)) then
            propName(1) = 'F1_N'
            propName(2) = 'F2_N'
            propName(3) = 'C_N'
            n(1:nbProp) = 20.d0
            unsurn(1:nbProp) = 1.d0
        end if
    end if

! - Get parameters
    call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                ' ', 'META_VISC', 0, ' ', [0.d0], &
                nbProp, propName, propVale, propCode, 2)
    if (present(n)) then
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'META_VISC', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        do iProp = 1, nbProp
            n(iProp) = propVale(iProp)
            unsurn(iProp) = 1.d0/n(iProp)
        end do
    end if

! - Name of parameters
    if (metaType .eq. META_STEEL) then
        if (present(c)) then
            propName(1) = 'F1_C'
            propName(2) = 'F2_C'
            propName(3) = 'F3_C'
            propName(4) = 'F4_C'
            propName(5) = 'C_C'
            c(1:nbProp) = 0.d0
        end if
    elseif (metaType .eq. META_ZIRC) then
        if (present(c)) then
            propName(1) = 'F1_C'
            propName(2) = 'F2_C'
            propName(3) = 'C_C'
            c(1:nbProp) = 0.d0
        end if
    end if

! - Get parameters
    if (present(c)) then
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'META_VISC', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        do iProp = 1, nbProp
            c(iProp) = propVale(iProp)
        end do
    end if

! - Name of parameters
    if (metaType .eq. META_STEEL) then
        if (present(m)) then
            propName(1) = 'F1_M'
            propName(2) = 'F2_M'
            propName(3) = 'F3_M'
            propName(4) = 'F4_M'
            propName(5) = 'C_M'
            m(1:nbProp) = 20.d0
        end if
    elseif (metaType .eq. META_ZIRC) then
        if (present(m)) then
            propName(1) = 'F1_M'
            propName(2) = 'F2_M'
            propName(3) = 'C_M'
            m(1:nbProp) = 20.d0
        end if
    end if

! - Get parameters
    if (present(m)) then
        call rcvalb(fami, kpg, ksp, poum, jvMaterCode, &
                    ' ', 'META_VISC', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 2)
        do iProp = 1, nbProp
            m(iProp) = propVale(iProp)
        end do
    end if
!
end subroutine
