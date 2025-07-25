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
!
!
#include "asterf_types.h"
!
interface
    subroutine zggevx(balanc, jobvl, jobvr, sense, n,&
                      a, lda, b, ldb, alpha,&
                      beta, vl, ldvl, vr, ldvr,&
                      ilo, ihi, lscale, rscale, abnrm,&
                      bbnrm, rconde, rcondv, work, lwork,&
                      rwork, iwork, bwork, info)
        blas_int, intent(in) :: ldvr
        blas_int, intent(in) :: ldvl
        blas_int, intent(in) :: ldb
        blas_int, intent(in) :: lda
        character(len=1), intent(in) :: balanc
        character(len=1), intent(in) :: jobvl
        character(len=1), intent(in) :: jobvr
        character(len=1), intent(in) :: sense
        blas_int, intent(in) :: n
        complex(kind=8), intent(inout) :: a(lda, *)
        complex(kind=8), intent(inout) :: b(ldb, *)
        complex(kind=8), intent(out) :: alpha(*)
        complex(kind=8), intent(out) :: beta(*)
        complex(kind=8), intent(out) :: vl(ldvl, *)
        complex(kind=8), intent(out) :: vr(ldvr, *)
        blas_int, intent(out) :: ilo
        blas_int, intent(out) :: ihi
        real(kind=8), intent(out) :: lscale(*)
        real(kind=8), intent(out) :: rscale(*)
        real(kind=8), intent(out) :: abnrm
        real(kind=8), intent(out) :: bbnrm
        real(kind=8), intent(out) :: rconde(*)
        real(kind=8), intent(out) :: rcondv(*)
        complex(kind=8), intent(out) :: work(*)
        blas_int, intent(in) :: lwork
        real(kind=8), intent(out) :: rwork(*)
        blas_int ,intent(out) :: iwork(*)
        logical(kind=4), intent(out) :: bwork(*)
        blas_int, intent(out) :: info
    end subroutine zggevx
end interface
