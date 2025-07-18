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
subroutine umatPrepareStrain(neps, epsm, deps, stran, dstran, &
                             dfgrd0, dfgrd1)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/r8inir.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
!
    integer(kind=8), intent(in) :: neps
    real(kind=8), intent(in) :: epsm(neps), deps(neps)
    real(kind=8), intent(out) :: stran(neps), dstran(neps)
    real(kind=8), intent(out) :: dfgrd0(3, 3), dfgrd1(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour (UMAT)
!
! Prepare stran, dstran, dfgrd0 and dfgrd1 variables
!
! --------------------------------------------------------------------------------------------------
!
! In  neps             : number of components of strains
! In  epsm             : mechanical strains at beginning of current step time
! In  deps             : increment of mechanical strains during current step time
! Out stran            : mechanical strains at beginning of current step time for UMAT
! Out dstran           : increment of mechanical strains during current step time for UMAT
! Out dfgrd0(3, 3)     : gradient of transformation at the beginning of the current step time
! Out dfgrd1(3, 3)     : gradient of transformation at the end of the current step time
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    if (neps .eq. 6) then
        b_n = to_blas_int(neps)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, deps, b_incx, dstran, b_incy)
        b_n = to_blas_int(neps)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsm, b_incx, stran, b_incy)
! TRAITEMENT DES COMPOSANTES 4,5,6 : DANS UMAT, GAMMAXY,XZ,YZ
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, rac2, dstran(4), b_incx)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, rac2, stran(4), b_incx)
! CAS DES GRANDES DEFORMATIONS : ON VEUT F- ET F+ -> non traité
        call r8inir(9, 0.d0, dfgrd0, 1)
        call r8inir(9, 0.d0, dfgrd1, 1)
    else
        ASSERT(.false.)
    end if
!
end subroutine
