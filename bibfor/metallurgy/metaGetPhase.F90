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
subroutine metaGetPhase(fami, poum, ipg, ispg, meta_type, &
                        nb_phase, phase_, zcold_, zhot_, tole_bound_)
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
    integer(kind=8), intent(in) :: ipg, ispg
    integer(kind=8), intent(in) :: meta_type
    integer(kind=8), intent(in) :: nb_phase
    real(kind=8), optional, intent(out) :: phase_(*)
    real(kind=8), optional, intent(out) :: zcold_
    real(kind=8), optional, intent(out) :: zhot_
    real(kind=8), optional, intent(in) :: tole_bound_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility - Metallurgy
!
! Get phase
!
! --------------------------------------------------------------------------------------------------
!
! In  fami         : Gauss family for integration point rule
! In  poum         : '-' or '+' for parameters evaluation (previous or current temperature)
! In  ipg          : current point gauss
! In  ispg         : current "sous-point" gauss
! In  meta_type    : type of metallurgy
! In  nb_phase    : total number of phase (cold and hot)
! Out phase       : phase
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
    integer(kind=8) :: i_phase_c, i_phase, iret, nb_phase_c
    real(kind=8) :: zcold, zhot, phase(5), tole_bound
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nb_phase .le. 5)
    phase(1:5) = 0.d0
    nb_phase_c = nb_phase-1
    if (present(tole_bound_)) then
        tole_bound = tole_bound_
    else
        tole_bound = r8prem()
    end if
!
! - Set cold phase
!
    do i_phase = 1, nb_phase
        if (meta_type .eq. META_STEEL) then
            call rcvarc('F', steel(i_phase), poum, fami, ipg, &
                        ispg, phase(i_phase), iret)
            if (iret .eq. 1) then
                phase(i_phase) = 0.d0
            end if
        elseif (meta_type .eq. META_ZIRC) then
            call rcvarc('F', zirc(i_phase), poum, fami, ipg, &
                        ispg, phase(i_phase), iret)
            if (iret .eq. 1) then
                phase(i_phase) = 0.d0
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
! - Sum of cold phases
!
    zcold = 0.d0
    if (meta_type .eq. META_STEEL) then
        call rcvarc('F', steel(6), poum, fami, ipg, &
                    ispg, zcold, iret)
        if (iret .eq. 1) then
            zcold = 0.d0
        end if
    elseif (meta_type .eq. META_ZIRC) then
        do i_phase_c = 1, nb_phase_c
            zcold = zcold+phase(i_phase_c)
        end do
        if (zcold .le. tole_bound) then
            zcold = 0.d0
        end if
        if (zcold .ge. 1.d0) then
            zcold = 1.d0
        end if
    end if
!
! - Set hot phase
!
    zhot = phase(nb_phase)
!
    if (present(phase_)) then
        phase_(1:nb_phase) = phase(1:nb_phase)
    end if
    if (present(zcold_)) then
        zcold_ = zcold
    end if
    if (present(zhot_)) then
        zhot_ = zhot
    end if
!
end subroutine
