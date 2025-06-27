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
    subroutine zggev(jobvl, jobvr, n, a, lda,&
                     b, ldb, alpha, beta, vl,&
                     ldvl, vr, ldvr, work, lwork,&
                     rwork, info)
        blas_int, intent(in) :: ldvr
        blas_int, intent(in) :: ldvl
        blas_int, intent(in) :: ldb
        blas_int, intent(in) :: lda
        character(len=1), intent(in) :: jobvl
        character(len=1), intent(in) :: jobvr
        blas_int, intent(in) :: n
        complex(kind=8), intent(inout) :: a(lda, *)
        complex(kind=8), intent(inout) :: b(ldb, *)
        complex(kind=8), intent(out) :: alpha(*)
        complex(kind=8), intent(out) :: beta(*)
        complex(kind=8), intent(out) :: vl(ldvl, *)
        complex(kind=8), intent(out) :: vr(ldvr, *)
        complex(kind=8), intent(out) :: work(*)
        blas_int, intent(in) :: lwork
        real(kind=8), intent(out) :: rwork(*)
        blas_int, intent(out) :: info
    end subroutine zggev
end interface
