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
module MetallurgyZirc_Compute_module
! ==================================================================================================
    use Metallurgy_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: metaInitZircGetPhases, metaInitCheckTransitionTime, metatInitZircSetField
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/utmess.h"
#include "asterfort/Metallurgy_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! metaInitZircGetPhases
!
! Get initial phases for zirc
!
! --------------------------------------------------------------------------------------------------
    subroutine metaInitZircGetPhases(jvPhaseIn, phase_tot)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: jvPhaseIn
        real(kind=8), intent(out) :: phase_tot
! ----- Local
        integer(kind=8) :: iPhase
!   ------------------------------------------------------------------------------------------------
!
        phase_tot = 0.d0
        do iPhase = 1, PZIRC_NB
            if (zr(jvPhaseIn-1+iPhase) .eq. r8vide() .or. isnan(zr(jvPhaseIn-1+iPhase))) then
                call utmess('F', 'META1_45')
            end if
            phase_tot = phase_tot+zr(jvPhaseIn-1+iPhase)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metaInitCheckTransitionTime
!
! Check transition time
!
! --------------------------------------------------------------------------------------------------
    subroutine metaInitCheckTransitionTime(nbPhase, jvPhaseIn)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbPhase, jvPhaseIn
!   ------------------------------------------------------------------------------------------------
!
        if (zr(jvPhaseIn-1+nbPhase+TIME_TRAN) .eq. r8vide() .or. &
            isnan(zr(jvPhaseIn-1+nbPhase+TIME_TRAN))) then
            call utmess('F', 'META1_47')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! metatInitZircSetField
!
! Set initial field for metallurgy
!
! --------------------------------------------------------------------------------------------------
    subroutine metatInitZircSetField(nbPhase, nbNode, nbVari, nbNodeMaxi, nbVariZirc, &
                                     jvTemp, jvPhaseIn, jvPhaseOut)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nbPhase, nbNode, nbVari
        integer(kind=8), intent(in) :: nbNodeMaxi, nbVariZirc
        integer(kind=8), intent(in) :: jvTemp, jvPhaseIn, jvPhaseOut
! ----- Local
        integer(kind=8) :: iNode, iVari
        real(kind=8) :: metaZirc(nbNodeMaxi*nbVariZirc), temp0
        real(kind=8) :: zalpha, zbeta
!   ------------------------------------------------------------------------------------------------
!
        metaZirc = 0.d0
        do iNode = 1, nbNode
            temp0 = zr(jvTemp+iNode-1)
            metaZirc(nbVari*(iNode-1)+PALPHA1) = zr(jvPhaseIn-1+PALPHA1)
            metaZirc(nbVari*(iNode-1)+PALPHA2) = zr(jvPhaseIn-1+PALPHA2)
            metaZirc(nbVari*(iNode-1)+nbPhase+ZIRC_TEMP) = temp0
            metaZirc(nbVari*(iNode-1)+nbPhase+TIME_TRAN) = zr(jvPhaseIn-1+nbPhase+TIME_TRAN)
            zalpha = metaZirc(nbVari*(iNode-1)+PALPHA1)+metaZirc(nbVari*(iNode-1)+PALPHA2)
            zbeta = 1.d0-zalpha
            if (zbeta .gt. 0.1d0) then
                metaZirc(nbVari*(iNode-1)+PALPHA1) = 0.d0
            else
                metaZirc(nbVari*(iNode-1)+PALPHA1) = 10.d0*(zalpha-0.9d0)*zalpha
            end if
            metaZirc(nbVari*(iNode-1)+PALPHA2) = zalpha- &
                                                 metaZirc(nbVari*(iNode-1)+PALPHA1)
            metaZirc(nbVari*(iNode-1)+PBETA) = zbeta
            do iVari = 1, nbVari
                zr(jvPhaseOut+nbVari*(iNode-1)-1+iVari) = metaZirc(nbVari*(iNode-1)+iVari)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module MetallurgyZirc_Compute_module
