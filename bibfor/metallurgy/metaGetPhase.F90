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
subroutine metaGetPhase(fami, poum, kpg, ksp, &
                        metaType, nbPhases, &
                        phase_, zcold_, zhot_, tole_bound_)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/rcvarc.h"
#include "asterfort/Metallurgy_type.h"
!
    character(len=*), intent(in) :: fami
    character(len=1), intent(in) :: poum
    integer(kind=8), intent(in) :: kpg, ksp
    integer(kind=8), intent(in) :: metaType, nbPhases
    real(kind=8), optional, intent(out) :: phase_(*)
    real(kind=8), optional, intent(out) :: zcold_
    real(kind=8), optional, intent(out) :: zhot_
    real(kind=8), optional, intent(in) :: tole_bound_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get phases
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  poum         : '-' or '+' for parameters evaluation (previous or current temperature)
! In  kpg          : current point gauss
! In  ksp          : current "sous-point" gauss
! In  metaType     : type of metallurgy
! In  nbPhases     : total number of phase (cold and hot)
! Out phase        : phases
! Out zcold        : sum of cold phase
! Out zhot         : hot phase
! In  tole_bound   : tolerance to project phase proportion on boundary
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: steel(6) = (/'PFERRITE', 'PPERLITE', &
                                                'PBAINITE', 'PMARTENS', &
                                                'PAUSTENI', 'PCOLDSUM'/)
    character(len=8), parameter :: zirc(3) = (/'ALPHPUR ', 'ALPHBETA', &
                                               'BETA    '/)
    integer(kind=8) :: iPhaseCold, iPhase, iret, nbPhasesCold
    real(kind=8) :: zcold, zhot, phase(5), tole_bound
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbPhases .le. 5)
    ASSERT(nbPhases .gt. 0)
    phase = 0.d0
    if (present(tole_bound_)) then
        tole_bound = tole_bound_
    else
        tole_bound = r8prem()
    end if
    nbPhasesCold = nbPhases-1

! - Set cold phase
    do iPhase = 1, nbPhases
        if (metaType .eq. META_STEEL) then
            call rcvarc('F', steel(iPhase), poum, fami, kpg, &
                        ksp, phase(iPhase), iret)
            if (iret .eq. 1) then
                phase(iPhase) = 0.d0
            end if
        elseif (metaType .eq. META_ZIRC) then
            call rcvarc('F', zirc(iPhase), poum, fami, kpg, &
                        ksp, phase(iPhase), iret)
            if (iret .eq. 1) then
                phase(iPhase) = 0.d0
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end do

! - Sum of cold phases
    zcold = 0.d0
    if (metaType .eq. META_STEEL) then
        call rcvarc('F', steel(6), poum, fami, kpg, &
                    ksp, zcold, iret)
        if (iret .eq. 1) then
            zcold = 0.d0
        end if

    elseif (metaType .eq. META_ZIRC) then
        do iPhaseCold = 1, nbPhasesCold
            zcold = zcold+phase(iPhaseCold)
        end do
        if (zcold .le. tole_bound) then
            zcold = 0.d0
        end if
        if (zcold .ge. 1.d0) then
            zcold = 1.d0
        end if
    end if

! - Set hot phase
    zhot = phase(nbPhases)
    if (present(phase_)) then
        phase_(1:nbPhases) = phase(1:nbPhases)
    end if
    if (present(zcold_)) then
        zcold_ = zcold
    end if
    if (present(zhot_)) then
        zhot_ = zhot
    end if
!
end subroutine
