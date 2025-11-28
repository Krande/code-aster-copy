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
subroutine verifepsa(famiZ, kpg, ksp, poumZ, epsiAnel)
!
    use BehaviourStrain_type
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/rcvarc.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: famiZ
    integer(kind=8), intent(in) :: kpg, ksp
    character(len=*), intent(in) :: poumZ
    real(kind=8), intent(out) :: epsiAnel(VARC_EPSA_NBCMP)
!
! --------------------------------------------------------------------------------------------------
!
! Get anelastic deformation (defined as external state variable)
!
! --------------------------------------------------------------------------------------------------
!
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! Out epsiAnel         : (generic) inelastic strains
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iCmp
    integer(kind=8) :: iretAnelPrev(VARC_EPSA_NBCMP), iretAnelCurr(VARC_EPSA_NBCMP)
    real(kind=8) :: epsiAnelPrev(VARC_EPSA_NBCMP), epsiAnelCurr(VARC_EPSA_NBCMP)
    character(len=1) :: poum
!
! --------------------------------------------------------------------------------------------------
!
    epsiAnel = 0.d0
    poum = poumZ

! - Get anelastic strains
    if (poum .eq. 'T' .or. poum .eq. '-') then
        epsiAnelPrev = 1
        epsiAnelPrev = r8nnem()
        do iCmp = 1, VARC_EPSA_NBCMP
            call rcvarc(' ', varcAnelName(iCmp), '-', famiZ, kpg, &
                        ksp, epsiAnelPrev(iCmp), iretAnelPrev(iCmp))
            if (iretAnelPrev(iCmp) .ne. 0) then
                epsiAnelPrev(iCmp) = 0.d0
            end if
        end do
    end if
    epsiAnelCurr = 1
    epsiAnelCurr = r8nnem()
    if (poum .eq. 'T' .or. poum .eq. '+') then
        do iCmp = 1, VARC_EPSA_NBCMP
            call rcvarc(' ', varcAnelName(iCmp), '+', famiZ, kpg, &
                        ksp, epsiAnelCurr(iCmp), iretAnelCurr(iCmp))
            if (iretAnelCurr(iCmp) .ne. 0) then
                epsiAnelCurr(iCmp) = 0.d0
            end if
        end do
    end if

! - Compute strains
    if (poum .eq. 'T') then
        do iCmp = 1, VARC_EPSA_NBCMP
            if (iretAnelPrev(iCmp)+iretAnelCurr(iCmp) .eq. 0) then
                epsiAnel(iCmp) = epsiAnelCurr(iCmp)-epsiAnelPrev(iCmp)
            end if
        end do
    else if (poum .eq. '-') then
        do iCmp = 1, VARC_EPSA_NBCMP
            if (iretAnelPrev(iCmp) .eq. 0) then
                epsiAnel(iCmp) = epsiAnelPrev(iCmp)
            end if
        end do
    else if (poum .eq. '+') then
        do iCmp = 1, VARC_EPSA_NBCMP
            if (iretAnelCurr(iCmp) .eq. 0) then
                epsiAnel(iCmp) = epsiAnelCurr(iCmp)
            end if
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
