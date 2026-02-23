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
subroutine te0320(option, nomte)
!
    use MetallurgySteel_Compute_module, only: metaInitSteelGetPhases, metaInitSteelCheckGrainSize
    use MetallurgySteel_Compute_module, only: metaInitSteelSumCold, metaInitSteelGetMartensite
    use MetallurgySteel_Compute_module, only: metatInitSteelSetField, metaSteelCheckFieldSize
    use MetallurgyZirc_Compute_module, only: metaInitCheckTransitionTime
    use MetallurgyZirc_Compute_module, only: metaInitZircGetPhases, metatInitZircSetField
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/Metallurgy_type.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: THERMIQUE - AXIS*, PLAN*
! Option: META_INIT_ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), pointer :: comporMeta(:) => null()
    character(len=16) :: metaName
    real(kind=8) :: ms0, phaseSum, phaseColdSum, grainSize
    integer(kind=8) :: nbNode, nbPhase, nbVari
    integer(kind=8) :: jvMater, jvTemp, jvPhaseIn, jvPhaseOut
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nbNode)
    ASSERT(nbNode .le. MT_NNOMAX2D)

! - Input fields
    call jevech('PMATERC', 'L', jvMater)
    call jevech('PTEMPER', 'L', jvTemp)
    call jevech('PCOMPME', 'L', vk16=comporMeta)
    call jevech('PPHASII', 'L', jvPhaseIn)

! - Output fields
    call jevech('PPHASOUT', 'E', jvPhaseOut)

! - Parameters from map
    metaName = comporMeta(ZMETATYPE)
    read (comporMeta(ZNBPHASE), '(I16)') nbPhase
    read (comporMeta(ZNBVARI), '(I16)') nbVari

! - Check size of field
    call metaSteelCheckFieldSize(metaName, jvPhaseIn)

! - Values required for META_INIT_ELNO vector
    phaseSum = 0.d0

    if (metaName .eq. 'ACIER') then
! ----- Get all phases
        call metaInitSteelGetPhases(metaName, jvPhaseIn, phaseSum)

! ----- Check size of grain
        call metaInitSteelCheckGrainSize(jvPhaseIn, grainSize)

! ----- Compute sum of "cold" phases
        call metaInitSteelSumCold(metaName, jvPhaseIn, phaseSum, phaseColdSum)

! ----- Get martensite temperature
        call metaInitSteelGetMartensite(jvMater, ms0)

    else if (metaName .eq. 'ZIRC') then
! ----- Get all phases
        call metaInitZircGetPhases(jvPhaseIn, phaseSum)

! ----- Check transition time
        call metaInitCheckTransitionTime(nbPhase, jvPhaseIn)

    else if (metaName .eq. 'VIDE') then

    else
        WRITE (6, *) "METANAME: ", metaName
        ASSERT(ASTER_FALSE)

    end if
!
    if (abs(phaseSum-1.d0) .gt. 1.d-2) then
        call utmess('F', 'META1_48', sr=phaseSum)
    end if

! - Set META_INI_ELNO field
    if (metaName .eq. 'ACIER') then
        ASSERT(nbVari .eq. NBVARISTEEL)
        call metatInitSteelSetField(nbNode, nbVari, MT_NNOMAX2D, NBVARISTEEL, &
                                    jvTemp, jvPhaseIn, jvPhaseOut, &
                                    ms0, phaseColdSum, grainSize)

    else if (metaName .eq. 'ZIRC') then
        ASSERT(nbVari .eq. NBVARIZIRC)
        call metatInitZircSetField(nbPhase, nbNode, nbVari, MT_NNOMAX2D, NBVARIZIRC, &
                                   jvTemp, jvPhaseIn, jvPhaseOut)

    else if (metaName .eq. 'VIDE') then

    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
