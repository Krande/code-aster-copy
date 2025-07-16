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
module MetallurgySteel_Compute_module
! ==================================================================================================
    use Metallurgy_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: metaInitSteelGetPhases, metaInitSteelCheckGrainSize, metaInitSteelSumCold
    public :: metaInitSteelGetMartensite, metatInitSteelSetField
    public :: metaSteelCheckFieldSize
    public :: metaSteelGetParameters, metaSteelTemperGetParameters, metaSteelTRCGetParameters
    public :: metaHardnessGetParameters
    public :: metaSteelCooling, metaSteelHeating, metaSteelCheckPhases
! ==================================================================================================
    private
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/metaSteelGrainSize.h"
#include "asterfort/rcadma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/smcarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
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
        integer(kind=8), intent(in) :: jvPhaseIn
        real(kind=8), intent(out) :: phase_tot
! ----- Local
        integer(kind=8) :: iPhase
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
        elseif (metaType .eq. 'ACIER_REVENU') then
            do iPhase = 1, PRSTEEL_NB
                if (zr(jvPhaseIn-1+iPhase) .eq. r8vide() .or. isnan(zr(jvPhaseIn-1+iPhase))) then
                    call utmess('F', 'META1_43')
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
        integer(kind=8), intent(in) :: nbPhase, jvPhaseIn
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
        integer(kind=8), intent(in) :: jvPhaseIn
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
        elseif (metaType .eq. 'ACIER_REVENU') then
            phase_scold = phase_tot-zr(jvPhaseIn-1+PRAUSTENITE)
            phase_ucold = zr(jvPhaseIn-1+PRSUMCOLD)
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
        integer(kind=8), intent(in) :: jvMater
        real(kind=8), intent(out) :: ms0
! ----- Local
        integer(kind=8), parameter :: kpg = 1, spt = 1
        character(len=8), parameter :: fami = 'FPG1', poum = '+'
        character(len=24), parameter :: paraName = "MS0"
        integer(kind=8) :: codret(1)
        real(kind=8) :: paraVale(1)
!   ------------------------------------------------------------------------------------------------
!
        call rcvalb(fami, kpg, spt, poum, zi(jvMater), &
                    ' ', 'META_ACIER', 0, ' ', [0.d0], &
                    1, paraName, paraVale, codret, 1)
        ms0 = paraVale(1)
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
        integer(kind=8), intent(in) :: nbPhase, nbNode, nbVari
        integer(kind=8), intent(in) :: nbNodeMaxi, nbVariSteel
        integer(kind=8), intent(in) :: jvTemp, jvPhaseIn, jvPhaseOut
        real(kind=8), intent(in) :: ms0, phase_ucold
! ----- Local
        integer(kind=8) :: iNode, iVari, iPhase
        real(kind=8) :: metaSteel(nbNodeMaxi*nbVariSteel), temp0
!   ------------------------------------------------------------------------------------------------
!
        metaSteel = 0.d0
        if (metaType .eq. 'ACIER') then
            do iNode = 1, nbNode
                temp0 = zr(jvTemp+iNode-1)
                do iPhase = 1, PSTEEL_NB
                    metaSteel(nbVari*(iNode-1)+iPhase) = zr(jvPhaseIn-1+iPhase)
                end do
                metaSteel(nbVari*(iNode-1)+PSUMCOLD) = phase_ucold
                metaSteel(nbVari*(iNode-1)+nbPhase+SIZE_GRAIN) = zr(jvPhaseIn-1+nbPhase+SIZE_GRAIN)
                metaSteel(nbVari*(iNode-1)+nbPhase+TEMP_MARTENSITE) = ms0
                metaSteel(nbVari*(iNode-1)+nbPhase+STEEL_TEMP) = temp0
                ASSERT(nbVari .eq. PSTEEL_NB+1+3)
                do iVari = 1, nbVari
                    zr(jvPhaseOut+nbVari*(iNode-1)-1+iVari) = metaSteel(nbVari*(iNode-1)+iVari)
                end do
            end do
        elseif (metaType .eq. 'ACIER_REVENU') then
            do iNode = 1, nbNode
                temp0 = zr(jvTemp+iNode-1)
                do iPhase = 1, PRSTEEL_NB
                    metaSteel(nbVari*(iNode-1)+iPhase) = zr(jvPhaseIn-1+iPhase)
                end do
                metaSteel(nbVari*(iNode-1)+PRSUMCOLD) = phase_ucold
                metaSteel(nbVari*(iNode-1)+nbPhase+SIZE_GRAIN) = zr(jvPhaseIn-1+nbPhase+SIZE_GRAIN)
                metaSteel(nbVari*(iNode-1)+nbPhase+TEMP_MARTENSITE) = ms0
                metaSteel(nbVari*(iNode-1)+nbPhase+STEEL_TEMP) = temp0
                metaSteel(nbVari*(iNode-1)+nbPhase+THER_CYCL) = 0.d0
                ASSERT(nbVari .eq. PRSTEEL_NB+1+4)
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
        integer(kind=8), intent(in) :: jvPhaseIn
! ----- Local
        integer(kind=8), parameter :: sizeFieldMaxi = 9
        integer(kind=8) :: iVari, fieldSize
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
        elseif (metaType .eq. 'ACIER_REVENU') then
            if (fieldSize .ne. PRVARIINIT) then
                call utmess("F", "META1_4", ni=2, vali=[PRVARIINIT, fieldSize])
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaSteelGetParameters
!
! Get parameters for steel behaviour in metallurgy
!
! --------------------------------------------------------------------------------------------------
    subroutine metaSteelGetParameters(jvMaterCode, metaSteelPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvMaterCode
        type(META_SteelParameters), intent(inout) :: metaSteelPara
! ----- Local
        integer(kind=8), parameter :: kpg = 1, spt = 1
        character(len=8), parameter :: fami = 'FPG1', poum = '+'
        integer(kind=8), parameter :: nbParaSteel = 7
        real(kind=8) :: paraSteelVale(nbParaSteel)
        integer(kind=8) :: codretSteel(nbParaSteel)
        character(len=16), parameter :: paraSteelName(nbParaSteel) = (/'AR3   ', &
                                                                       'ALPHA ', &
                                                                       'MS0   ', &
                                                                       'AC1   ', &
                                                                       'AC3   ', &
                                                                       'TAUX_1', &
                                                                       'TAUX_3'/)
        integer(kind=8), parameter :: nbParaAuste = 4
        real(kind=8) :: paraAusteVale(nbParaAuste)
        integer(kind=8) :: codretAuste(nbParaAuste)
        character(len=16), parameter :: paraAusteName(nbParaAuste) = (/'LAMBDA0', &
                                                                       'QSR_K  ', &
                                                                       'D10    ', &
                                                                       'WSR_K  '/)
!   ------------------------------------------------------------------------------------------------
!
        paraSteelVale = 0.d0
        call rcvalb(fami, kpg, spt, poum, &
                    jvMaterCode, ' ', 'META_ACIER', &
                    1, 'INST', [0.d0], &
                    nbParaSteel, paraSteelName, paraSteelVale, &
                    codretSteel, iarret=1)
        metaSteelPara%ar3 = paraSteelVale(1)
        metaSteelPara%alpha = paraSteelVale(2)
        metaSteelPara%ms0 = paraSteelVale(3)
        metaSteelPara%ac1 = paraSteelVale(4)
        metaSteelPara%ac3 = paraSteelVale(5)
        metaSteelPara%taux_1 = paraSteelVale(6)
        metaSteelPara%taux_3 = paraSteelVale(7)

        paraAusteVale = 0.d0
        call rcvalb(fami, kpg, spt, poum, &
                    jvMaterCode, ' ', 'META_ACIER', &
                    1, 'INST', [0.d0], &
                    nbParaAuste, paraAusteName, paraAusteVale, &
                    codretAuste, iarret=0, nan='NON')
        metaSteelPara%austenite%lambda0 = paraAusteVale(1)
        metaSteelPara%austenite%qsr_k = paraAusteVale(2)
        metaSteelPara%austenite%d10 = paraAusteVale(3)
        metaSteelPara%austenite%wsr_k = paraAusteVale(4)
        if ((codretAuste(1) .eq. 0) .and. (codretAuste(3) .eq. 1)) then
            call utmess('F', 'METALLURGY1_73')
        end if

! ----- Update size of martensite grain ?
        if (codretAuste(1) .eq. 0) then
            metaSteelPara%l_grain_size = ASTER_TRUE
            if (metaSteelPara%austenite%lambda0 .le. r8prem()) then
                metaSteelPara%l_grain_size = ASTER_FALSE
            end if
        else
            metaSteelPara%l_grain_size = ASTER_FALSE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaSteelTemperGetParameters
!
! Get parameters for steel with tempering behaviour in metallurgy
!
! --------------------------------------------------------------------------------------------------
    subroutine metaSteelTemperGetParameters(jvMaterCode, metaSteelPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvMaterCode
        type(META_SteelParameters), intent(inout) :: metaSteelPara
! ----- Local
        integer(kind=8), parameter :: kpg = 1, spt = 1
        character(len=8), parameter :: fami = 'FPG1', poum = '+'
        integer(kind=8), parameter :: nbParaTemper = 6
        real(kind=8) :: paraTemperVale(nbParaTemper)
        integer(kind=8) :: codretTemper(nbParaTemper)
        character(len=16), parameter :: paraTemperName(nbParaTemper) = (/'BAINITE_B    ', &
                                                                         'BAINITE_N    ', &
                                                                         'MARTENSITE_B ', &
                                                                         'MARTENSITE_N ', &
                                                                         'TEMP         ', &
                                                                         'TEMP_MAINTIEN'/)
!   ------------------------------------------------------------------------------------------------
!
        paraTemperVale = 0.d0
        call rcvalb(fami, kpg, spt, poum, &
                    jvMaterCode, ' ', 'META_ACIER_REVENU', &
                    1, 'INST', [0.d0], &
                    nbParaTemper, paraTemperName, paraTemperVale, &
                    codretTemper, iarret=0)
        if ((codretTemper(1) .eq. 0) .and. (codretTemper(2) .eq. 0) .and. &
            (codretTemper(3) .eq. 0) .and. (codretTemper(4) .eq. 0) .and. &
            (codretTemper(5) .eq. 0) .and. (codretTemper(6) .eq. 0)) then
            metaSteelPara%temper%bainite_b = paraTemperVale(1)
            metaSteelPara%temper%bainite_n = paraTemperVale(2)
            metaSteelPara%temper%martensite_b = paraTemperVale(3)
            metaSteelPara%temper%martensite_n = paraTemperVale(4)
            metaSteelPara%temper%temp = paraTemperVale(5)
            metaSteelPara%temper%tempHold = paraTemperVale(6)
        else
            call utmess('F', 'METALLURGY1_74')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaSteelTRCGetParameters
!
! Get parameters for TRC curves in metallurgy
!
! --------------------------------------------------------------------------------------------------
    subroutine metaSteelTRCGetParameters(jvMaterCode, metaSteelPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvMaterCode
        type(META_SteelParameters), intent(inout) :: metaSteelPara
! ----- Local
        real(kind=8), parameter :: toleTemp = 10.
        integer(kind=8) :: jvPftrc
        integer(kind=8) :: nbcb1, nbcb2, nblexp
        integer(kind=8) :: codret, nbTrc, nbHist, nbExp, iHist, iExp
        integer(kind=8) :: jftrc, jtrc
        integer(kind=8) :: iadexp, iadckm, iadtrc, shift
        real(kind=8) :: tempAR3FromMate, tempAR3FromTRC
        real(kind=8) :: tempPrev, tempCurr
        aster_logical :: lCooling
!   ------------------------------------------------------------------------------------------------
!
        call jevech('PFTRC', 'L', jvPftrc)
        jftrc = zi(jvPftrc)
        jtrc = zi(jvPftrc+1)
        call rcadma(jvMaterCode, 'META_ACIER', 'TRC', iadtrc, codret, 1)
        nbcb1 = nint(zr(iadtrc+1))
        nbHist = nint(zr(iadtrc+2))
        nbcb2 = nint(zr(iadtrc+1+2+nbcb1*nbHist))
        nblexp = nint(zr(iadtrc+1+2+nbcb1*nbHist+1))
        nbTrc = nint(zr(iadtrc+1+2+nbcb1*nbHist+2+nbcb2*nblexp+1))
        ASSERT(nbTrc .eq. 1)
        iadexp = 5+nbcb1*nbHist
        iadckm = 7+nbcb1*nbHist+nbcb2*nblexp
        metaSteelPara%trc%jv_ftrc = jftrc
        metaSteelPara%trc%jv_trc = jtrc
        metaSteelPara%trc%iadexp = iadexp
        metaSteelPara%trc%iadtrc = iadtrc
        metaSteelPara%trc%nbHist = nbHist

! ----- Parameters for martensite evolution
        metaSteelPara%trc%martensiteLaw%austeniteMin = zr(iadtrc+iadckm-1+1)
        metaSteelPara%trc%martensiteLaw%akm = zr(iadtrc+iadckm-1+2)
        metaSteelPara%trc%martensiteLaw%bkm = zr(iadtrc+iadckm-1+3)
        metaSteelPara%trc%martensiteLaw%lowerSpeed = zr(iadtrc+iadckm-1+4)

! ----- Parameters for size of austenite grain
        metaSteelPara%trc%austeniteGrain%dref = zr(iadtrc+iadckm-1+5)
        metaSteelPara%trc%austeniteGrain%a = zr(iadtrc+iadckm-1+6)

! ----- Check consistency of temperature
        tempAR3FromMate = metaSteelPara%ar3
        shift = 0
        do iHist = 1, nbHist
            nbExp = nint(zr(iadtrc+11+9*(iHist-1)))
            lCooling = ASTER_FALSE
            do iExp = 1, nbExp-1
                tempPrev = zr(iadtrc+iadexp-1+4*(shift+iExp))
                tempCurr = zr(iadtrc+iadexp-1+4*(shift+iExp+1))
                if (iExp .ge. 2) then
                    if (tempCurr .le. tempPrev) then
                        lCooling = ASTER_TRUE
                        exit
                    end if
                end if
            end do
            if (lCooling) then
                iExp = 1
                tempAR3FromTRC = zr(iadtrc+iadexp-1+4*(shift+iExp))
                if (abs(tempAR3FromMate-tempAR3FromTRC) .gt. toleTemp) then
                    call utmess('A', "META1_50", sr=toleTemp)
                end if
            end if
            shift = shift+nbExp
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaHardnessGetParameters
!
! Get parameters for hardness un metallurgy
!
! --------------------------------------------------------------------------------------------------
    subroutine metaHardnessGetParameters(jvMaterCode, metaType, metaHardnessPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvMaterCode
        character(len=16), intent(in) :: metaType
        type(META_HardnessParameters), intent(inout) :: metaHardnessPara
! ----- Local
        integer(kind=8), parameter :: kpg = 1, spt = 1
        character(len=8), parameter :: fami = 'FPG1', poum = '+'
        integer(kind=8), parameter :: nbParaHard = PSTEEL_NB
        real(kind=8) :: paraHardVale(nbParaHard)
        integer(kind=8) :: codretHard(nbParaHard)
        character(len=16), parameter :: paraHardName(nbParaHard) = (/'F1_DURT', &
                                                                     'F2_DURT', &
                                                                     'F3_DURT', &
                                                                     'F4_DURT', &
                                                                     'C_DURT '/)
        integer(kind=8), parameter :: nbParaHardTemp = PRSTEEL_NB
        real(kind=8) :: paraHardTempVale(nbParaHardTemp)
        integer(kind=8) :: codretHardTemp(nbParaHardTemp)
        character(len=16), parameter :: paraHardTempName(nbParaHardTemp) = (/'F1_DURT       ', &
                                                                             'F2_DURT       ', &
                                                                             'F3_DURT       ', &
                                                                             'F3_REVENU_DURT', &
                                                                             'F4_DURT       ', &
                                                                             'F4_REVENU_DURT', &
                                                                             'C_DURT        '/)
!   ------------------------------------------------------------------------------------------------
!
        if (metaType .eq. "ACIER") then
            call rcvalb(fami, kpg, spt, poum, jvMaterCode, &
                        ' ', 'DURT_META', 1, 'TEMP', [0.d0], &
                        nbParaHard, paraHardName, paraHardVale, codretHard, 2)
            metaHardnessPara%hardSteel(1:PSTEEL_NB) = paraHardVale
        else if (metaType .eq. "ACIER_REVENU") then
            call rcvalb(fami, kpg, spt, poum, jvMaterCode, &
                        ' ', 'DURT_META', 1, 'TEMP', [0.d0], &
                        nbParaHardTemp, paraHardTempName, paraHardTempVale, codretHardTemp, 2)
            metaHardnessPara%hardSteel(1:PRSTEEL_NB) = paraHardTempVale
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaSteelCooling
!
! Effect of cooling on steel phases (austenite => ferrite)
!
! --------------------------------------------------------------------------------------------------
    subroutine metaSteelCooling(metaSteelPara, &
                                temp0, temp1, temp2, &
                                nbStep, deltaTime01, deltaTime12, &
                                nbVari, nbPhase, nbHist, &
                                metaPrev, metaCurr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbStep
        real(kind=8), intent(in) :: temp0, temp1, temp2
        real(kind=8), intent(in) :: deltaTime01, deltaTime12
        integer(kind=8), intent(in) :: nbVari, nbPhase, nbHist
        type(META_SteelParameters), intent(in) :: metaSteelPara
        real(kind=8), intent(in) :: metaPrev(nbVari)
        real(kind=8), intent(out) :: metaCurr(nbVari)
! ----- Local
        real(kind=8) :: tempStep, metaDumm(9), deltaTimeStep, tPoint
        integer(kind=8) :: iStep
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nbVari .eq. 9)
        metaDumm = metaPrev
        deltaTimeStep = deltaTime12/dble(nbStep)

        do iStep = 1, nbStep
            tempStep = temp1+(temp2-temp1)*dble(iStep-1)/dble(nbStep)
            metaCurr(nbPhase+STEEL_TEMP) = temp1+ &
                                           (dble(iStep)*(temp2-temp1))/dble(nbStep)
            if (nbStep .eq. 1) then
                tPoint = (temp1-temp0)/deltaTime01
            else
                tPoint = (metaCurr(nbPhase+STEEL_TEMP)-tempStep)/deltaTimeStep
            end if
            call smcarc(nbHist, nbPhase, &
                        zr(metaSteelPara%trc%jv_ftrc), &
                        zr(metaSteelPara%trc%jv_trc), &
                        zr(metaSteelPara%trc%iadtrc+3), &
                        zr(metaSteelPara%trc%iadtrc+metaSteelPara%trc%iadexp), &
                        metaSteelPara, &
                        tempStep, tPoint, deltaTime01, &
                        metaDumm, metaCurr)
            metaDumm(:) = metaCurr(:)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaSteelHeating
!
! Effect of heating on steel phases (austenite <= ferrite)
!
! --------------------------------------------------------------------------------------------------
    subroutine metaSteelHeating(metaSteelPara, &
                                temp0, temp1, temp2, &
                                nbStep, deltaTime01, deltaTime12, &
                                zAustPrev, &
                                dPrev, dCurr, zAustCurr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbStep
        real(kind=8), intent(in) :: temp0, temp1, temp2
        real(kind=8), intent(in) :: deltaTime01, deltaTime12
        real(kind=8), intent(inout) :: zAustPrev
        type(META_SteelParameters), intent(in) :: metaSteelPara
        real(kind=8), intent(in) :: dPrev
        real(kind=8), intent(out) :: dCurr, zAustCurr
! ----- Local
        real(kind=8) :: ti1, ti2, deltaTimeStep, tPoint
        real(kind=8) :: zeq1, zeq2, taux, coef_phase, dStep
        integer(kind=8) :: iStep
        real(kind=8), parameter :: un = 1.d0, zero = 0.d0
        real(kind=8), parameter :: epsi = 1.d-10
!   ------------------------------------------------------------------------------------------------
!
        deltaTimeStep = deltaTime12/dble(nbStep)

        dStep = dPrev
        do iStep = 1, nbStep
            if (nbStep .eq. 1) then
                ti1 = temp1
                ti2 = temp2
                tPoint = (temp1-temp0)/deltaTime01
            else
                ti1 = temp1+(temp2-temp1)*dble(iStep-1)/dble(nbStep)
                ti2 = temp1+(temp2-temp1)*dble(iStep)/dble(nbStep)
                tpoint = (ti2-ti1)/deltaTimeStep
            end if
            zeq1 = min((ti1-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
            zeq2 = min((ti2-metaSteelPara%ac1)/(metaSteelPara%ac3-metaSteelPara%ac1), un)
            taux = metaSteelPara%taux_1+(metaSteelPara%taux_3-metaSteelPara%taux_1)*zeq1

            if ((ti1 .lt. (metaSteelPara%ac1-epsi)) .or. (zAustPrev .ge. un)) then
                zAustCurr = zAustPrev
            else
                if (zeq2 .ge. (un-epsi)) then
                    tpoint = zero
                end if
                if (zAustPrev .gt. zeq1) then
                    zAustCurr = (taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))
                    zAustCurr = zAustCurr*exp(-deltaTimeStep/taux)
                    zAustCurr = ((-taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))+ &
                                 zeq2+zAustCurr-zeq1)* &
                                (un-zAustPrev)/(un-zeq1)
                    zAustCurr = zAustCurr+zAustPrev
                else
                    zAustCurr = (taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))-zeq1+zAustPrev
                    zAustCurr = zAustCurr*exp(-deltaTimeStep/taux)
                    zAustCurr = (-taux*tpoint/(metaSteelPara%ac3-metaSteelPara%ac1))+zeq2+zAustCurr
                end if
            end if

! --------- Compute size of grain
            coef_phase = 1.d0
            if (abs(zAustCurr) .ge. epsi) then
                coef_phase = (1.d0-(zAustCurr-zAustPrev)/zAustCurr)
            end if
            call metaSteelGrainSize(metaSteelPara, &
                                    ti1, deltaTime01, deltaTime12, &
                                    zAustCurr, coef_phase, &
                                    dStep, dCurr)
            if (metaSteelPara%l_grain_size .and. nbStep .gt. 1) then
                zAustPrev = zAustCurr
                dStep = dCurr
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaSteelCheckPhases
!
! Check domain of definition of phases
!
! --------------------------------------------------------------------------------------------------
    subroutine metaSteelCheckPhases(nbVari, metaCurr)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbVari
        real(kind=8), intent(in) :: metaCurr(nbVari)
! ----- Local
        integer(kind=8) :: iPhase
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nbVari .ge. PSTEEL_NB)
        do iPhase = 1, PSTEEL_NB
            if (metaCurr(iPhase) .lt. 0.d0 .or. metaCurr(iPhase) .gt. 1.d0) then
                call utmess('F', 'META1_52', si=iPhase, sr=metaCurr(iPhase))
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module MetallurgySteel_Compute_module
