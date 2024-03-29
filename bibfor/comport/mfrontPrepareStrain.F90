! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine mfrontPrepareStrain(l_greenlag, l_pred, &
                               neps, epsm, deps, &
                               stran, dstran)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcdetf.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
!
    aster_logical, intent(in) :: l_greenlag, l_pred
    integer, intent(in) :: neps
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
! In  l_pred           : flag if prediction
! In  option           : option of calcul : RIGI_MECA, FULL_MECA...
! In  neps             : number of components of strains
! In  epsm             : mechanical strains at T- for all kinematics but simo_miehe
!                        total strains at T- for simo_miehe
! In  deps             : incr of mechanical strains during step time for all but simo_miehe
!                        incr of total strains during step time for simo_miehe
! Out stran            : mechanical strains at beginning of current step time for MFront
! Out dstran           : increment of mechanical strains during step time for MFront
! Out detf             : determinant of gradient
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: rac2 = 1.0
    !sqrt(2.d0)
    real(kind=8) :: dfgrd0(3, 3), dfgrd1(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
    stran(1:neps) = 0.d0
    dstran(1:neps) = 0.d0
!
    if (l_greenlag) then
        ASSERT(neps .eq. 9)
        dfgrd0(:, :) = 0.d0
        dfgrd1(:, :) = 0.d0
        call dcopy(neps, epsm, 1, dfgrd0, 1)
        if (l_pred) then
            call dcopy(neps, epsm, 1, dfgrd1, 1)
        else
            call dcopy(neps, deps, 1, dfgrd1, 1)
        end if
        call dcopy(neps, dfgrd0, 1, stran, 1)
        call dcopy(neps, dfgrd1, 1, dstran, 1)
    else
        ASSERT(neps .ne. 9)
        call dcopy(neps, deps, 1, dstran, 1)
        call dcopy(neps, epsm, 1, stran, 1)
        if (neps .eq. 6) then
            call dscal(3, rac2, dstran(4), 1)
            call dscal(3, rac2, stran(4), 1)
        end if
        if (neps .eq. 4) then
            call dscal(1, rac2, dstran(4), 1)
            call dscal(1, rac2, stran(4), 1)
        end if
    end if
!
end subroutine
