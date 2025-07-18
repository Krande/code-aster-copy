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
    subroutine zlatrs(uplo, trans, diag, normin, n,&
                      a, lda, x, scale, cnorm,&
                      info)
        blas_int, intent(in) :: lda
        character(len=1), intent(in) :: uplo
        character(len=1), intent(in) :: trans
        character(len=1), intent(in) :: diag
        character(len=1), intent(in) :: normin
        blas_int, intent(in) :: n
        complex(kind=8), intent(in) :: a(lda, *)
        complex(kind=8), intent(inout) :: x(*)
        real(kind=8), intent(out) :: scale
        real(kind=8), intent(inout) :: cnorm(*)
        blas_int, intent(out) :: info
    end subroutine zlatrs
end interface
