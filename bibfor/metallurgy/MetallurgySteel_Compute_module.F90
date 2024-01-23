! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
module MetallurgySteel_Compute_module
! ==================================================================================================
    use Metallurgy_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: metaInitSteelGetPhases, metaInitSteelCheckGrainSize, metaInitSteelSumCold
    public :: metaInitSteelGetMartensite, metatInitSteelSetField
    public :: metaSteelCheckFieldSize
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/Metallurgy_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! metaInitSteelGetPhases
!
! Get initial phases for steel
!
! --------------------------------------------------------------------------------------------------
    subroutine metaInitSteelGetPhases(metaType, jvPhaseIn, phase_tot)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: metaType
        integer, intent(in) :: jvPhaseIn
        real(kind=8), intent(out) :: phase_tot
! ----- Local
        integer :: iPhase
!   ------------------------------------------------------------------------------------------------
!
        phase_tot = 0.d0
        if (metaType .eq. 'ACIER') then
            do iPhase = 1, PSTEEL_NB
                if (zr(jvPhaseIn-1+iPhase) .eq. r8vide() .or. isnan(zr(jvPhaseIn-1+iPhase))) then
                    call utmess('F', 'META1_44')
                end if
                phase_tot = phase_tot+zr(jvPhaseIn-1+iPhase)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaInitSteelCheckGrainSize
!
! Check size of grain
!
! --------------------------------------------------------------------------------------------------
    subroutine metaInitSteelCheckGrainSize(nbPhase, jvPhaseIn)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: nbPhase, jvPhaseIn
!   ------------------------------------------------------------------------------------------------
!
        if (zr(jvPhaseIn-1+nbPhase+SIZE_GRAIN) .eq. r8vide() .or. &
            isnan(zr(jvPhaseIn-1+nbPhase+SIZE_GRAIN))) then
            call utmess('F', 'META1_46')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaInitSteelSumCold
!
! Compute sum of cold phases
!
! --------------------------------------------------------------------------------------------------
    subroutine metaInitSteelSumCold(metaType, jvPhaseIn, phase_tot, phase_ucold)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: metaType
        integer, intent(in) :: jvPhaseIn
        real(kind=8), intent(in) :: phase_tot
        real(kind=8), intent(out) :: phase_ucold
! ----- Local
        real(kind=8) :: phase_scold
!   ------------------------------------------------------------------------------------------------
!
        phase_ucold = 0.d0
        if (metaType .eq. 'ACIER') then
            phase_scold = phase_tot-zr(jvPhaseIn-1+PAUSTENITE)
            phase_ucold = zr(jvPhaseIn-1+PSUMCOLD)
            if (abs(phase_scold-phase_ucold) .gt. 1.d-2 .or. phase_ucold .eq. r8vide()) then
                call utmess('A', 'META1_49')
                phase_ucold = phase_scold
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaInitSteelGetMartensite
!
! Get martensite temperature
!
! --------------------------------------------------------------------------------------------------
    subroutine metaInitSteelGetMartensite(jvMater, ms0)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer, intent(in) :: jvMater
        real(kind=8), intent(out) :: ms0
! ----- Local
        integer, parameter :: kpg = 1, spt = 1
        character(len=8), parameter :: fami = 'FPG1', poum = '+'
        character(len=24), parameter :: resuName = "MS0"
        integer :: icodre(1)
        real(kind=8) :: resuVale(1)
!   ------------------------------------------------------------------------------------------------
!
        call rcvalb(fami, kpg, spt, poum, zi(jvMater), &
                    ' ', 'META_ACIER', 0, ' ', [0.d0], &
                    1, resuName, resuVale, icodre, 1)
        ms0 = resuVale(1)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metatInitSteelSetField
!
! Set initial field for metallurgy
!
! --------------------------------------------------------------------------------------------------
    subroutine metatInitSteelSetField(metaType, &
                                      nbPhase, nbNode, nbVari, nbNodeMaxi, nbVariSteel, &
                                      jvTemp, jvPhaseIn, jvPhaseOut, &
                                      ms0, phase_ucold)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: metaType
        integer, intent(in) :: nbPhase, nbNode, nbVari
        integer, intent(in) :: nbNodeMaxi, nbVariSteel
        integer, intent(in) :: jvTemp, jvPhaseIn, jvPhaseOut
        real(kind=8), intent(in) :: ms0, phase_ucold
! ----- Local
        integer :: iNode, iVari, iPhase
        real(kind=8) :: metaSteel(nbNodeMaxi*nbVariSteel), tno0
!   ------------------------------------------------------------------------------------------------
!
        metaSteel = 0.d0
        if (metaType .eq. 'ACIER') then
            do iNode = 1, nbNode
                tno0 = zr(jvTemp+iNode-1)
                do iPhase = 1, PSTEEL_NB
                    metaSteel(nbVari*(iNode-1)+iPhase) = zr(jvPhaseIn-1+iPhase)
                end do
                metaSteel(nbVari*(iNode-1)+PSUMCOLD) = phase_ucold
                metaSteel(nbVari*(iNode-1)+nbPhase+SIZE_GRAIN) = zr(jvPhaseIn-1+nbPhase+SIZE_GRAIN)
                metaSteel(nbVari*(iNode-1)+nbPhase+TEMP_MARTENSITE) = ms0
                metaSteel(nbVari*(iNode-1)+nbPhase+STEEL_TEMP) = tno0
                ASSERT(nbVari .eq. PSTEEL_NB+1+3)
                do iVari = 1, nbVari
                    zr(jvPhaseOut+nbVari*(iNode-1)-1+iVari) = metaSteel(nbVari*(iNode-1)+iVari)
                end do
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaSteelCheckFieldSize
!
! Check size of field
!
! --------------------------------------------------------------------------------------------------
    subroutine metaSteelCheckFieldSize(metaType, jvPhaseIn)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: metaType
        integer, intent(in) :: jvPhaseIn
! ----- Local
        integer, parameter :: sizeFieldMaxi = 9
        integer :: iVari, fieldSize
!   ------------------------------------------------------------------------------------------------
!
        fieldSize = 0
        do iVari = 1, sizeFieldMaxi
            if (zr(jvPhaseIn-1+iVari) .ne. r8vide()) then
                fieldSize = fieldSize+1
            end if
        end do
        if (metaType .eq. 'ACIER') then
            if (fieldSize .lt. PVARIINIT) then
                call utmess("F", "META1_4", ni=2, vali=[PVARIINIT, fieldSize])
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module MetallurgySteel_Compute_module
