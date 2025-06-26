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
subroutine mfrontPrepareStrain(l_greenlag, neps, epsm, deps, stran, dstran)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
!
    aster_logical, intent(in) :: l_greenlag
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps), deps(neps)
    real(kind=8), intent(out) :: stran(neps), dstran(neps)
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour (MFront)
!
! Prepare transformation gradient for large strains
! Prepare stran and dstran
!
! --------------------------------------------------------------------------------------------------
!
! In  l_greenlag       : .true. if large strains with GREEN_LAGRANGE
! In  neps             : number of components of strains
! In  epsm             : mechanical strains at T- for all kinematics but simo_miehe
!                        total strains at T- for simo_miehe
! In  deps             : incr of mechanical strains during step time for all but simo_miehe
!                        incr of total strains during step time for simo_miehe
! Out stran            : mechanical strains at beginning of current step time for MFront
! Out dstran           : increment of mechanical strains during step time for MFront
!
! --------------------------------------------------------------------------------------------------
!
    stran = 0.d0
    dstran = 0.d0
!
    if (l_greenlag) then
        ASSERT(neps .eq. 9)
        ! - Reordering to be consistent with mfront
        stran(1) = epsm(1)
        stran(2) = epsm(5)
        stran(3) = epsm(9)
        stran(4) = epsm(4)
        stran(5) = epsm(2)
        stran(6) = epsm(7)
        stran(7) = epsm(3)
        stran(8) = epsm(8)
        stran(9) = epsm(6)
        ! - Reordering to be consistent with mfront
        dstran(1) = deps(1)
        dstran(2) = deps(5)
        dstran(3) = deps(9)
        dstran(4) = deps(4)
        dstran(5) = deps(2)
        dstran(6) = deps(7)
        dstran(7) = deps(3)
        dstran(8) = deps(8)
        dstran(9) = deps(6)
    else
        ASSERT(neps .ne. 9)
        stran = epsm
        dstran = deps
    end if
!
end subroutine
