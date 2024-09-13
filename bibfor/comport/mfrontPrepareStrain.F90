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
subroutine mfrontPrepareStrain(l_greenlag, l_pred, neps, epsm, deps, &
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
    blas_int :: b_incx, b_incy, b_n
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
        b_n = to_blas_int(neps)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsm, b_incx, dfgrd0, b_incy)
        if (l_pred) then
            b_n = to_blas_int(neps)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, epsm, b_incx, dfgrd1, b_incy)
        else
            b_n = to_blas_int(neps)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, deps, b_incx, dfgrd1, b_incy)
        end if
        b_n = to_blas_int(neps)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dfgrd0, b_incx, stran, b_incy)
        b_n = to_blas_int(neps)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dfgrd1, b_incx, dstran, b_incy)
    else
        ASSERT(neps .ne. 9)
        b_n = to_blas_int(neps)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, deps, b_incx, dstran, b_incy)
        b_n = to_blas_int(neps)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsm, b_incx, stran, b_incy)
        if (neps .eq. 6) then
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, rac2, dstran(4), b_incx)
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            call dscal(b_n, rac2, stran(4), b_incx)
        end if
        if (neps .eq. 4) then
            b_n = to_blas_int(1)
            b_incx = to_blas_int(1)
            call dscal(b_n, rac2, dstran(4), b_incx)
            b_n = to_blas_int(1)
            b_incx = to_blas_int(1)
            call dscal(b_n, rac2, stran(4), b_incx)
        end if
    end if
!
end subroutine
