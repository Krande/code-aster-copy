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
! aslint: disable=W0413
!
subroutine dpassa(anglNaut, irep, matr_tran)
!
    implicit none
!
#include "asterfort/matrot.h"
!
    real(kind=8), intent(in) :: anglNaut(3)
    integer(kind=8), intent(out) :: irep
    real(kind=8), intent(out) :: matr_tran(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
! Elasticity
!
! Construct transition matrix for orthotropic elasticity - 3D case
!
! --------------------------------------------------------------------------------------------------
!
! In  anglNaut         : nautical angles
! Out irep             : 0 if matrix is trivial (identity), 1 otherwise
! Out matr_tran        : transition matrix
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: p(3, 3)
    real(kind=8), parameter :: deux = 2.d0, zero = 0.d0
!
! --------------------------------------------------------------------------------------------------
!
    irep = 0
    p = 0.d0

! - Compute matrix from orthotropic reference frame to global reference frame
    if (anglNaut(1) .eq. zero .and. anglNaut(2) .eq. zero .and. anglNaut(3) .eq. zero) then
        irep = 0
    else
        call matrot(anglNaut, p)
        irep = 1
    end if

! - Transition matrix for elasticity (fourth order)
! ---- CETTE MATRICE EST CONSTRUITE EN PARTANT DE LA CONSIDERATION QUE
! ----  (SIGMA_GLOB):(EPSILON_GLOB) = (SIGMA_ORTH):(EPSILON_ORTH)
    if (irep .eq. 1) then
        matr_tran(1, 1) = p(1, 1)*p(1, 1)
        matr_tran(1, 2) = p(1, 2)*p(1, 2)
        matr_tran(1, 3) = p(1, 3)*p(1, 3)
        matr_tran(1, 4) = p(1, 1)*p(1, 2)
        matr_tran(1, 5) = p(1, 1)*p(1, 3)
        matr_tran(1, 6) = p(1, 2)*p(1, 3)
!
        matr_tran(2, 1) = p(2, 1)*p(2, 1)
        matr_tran(2, 2) = p(2, 2)*p(2, 2)
        matr_tran(2, 3) = p(2, 3)*p(2, 3)
        matr_tran(2, 4) = p(2, 1)*p(2, 2)
        matr_tran(2, 5) = p(2, 1)*p(2, 3)
        matr_tran(2, 6) = p(2, 2)*p(2, 3)
!
        matr_tran(3, 1) = p(3, 1)*p(3, 1)
        matr_tran(3, 2) = p(3, 2)*p(3, 2)
        matr_tran(3, 3) = p(3, 3)*p(3, 3)
        matr_tran(3, 4) = p(3, 1)*p(3, 2)
        matr_tran(3, 5) = p(3, 1)*p(3, 3)
        matr_tran(3, 6) = p(3, 2)*p(3, 3)
!
        matr_tran(4, 1) = deux*p(1, 1)*p(2, 1)
        matr_tran(4, 2) = deux*p(1, 2)*p(2, 2)
        matr_tran(4, 3) = deux*p(1, 3)*p(2, 3)
        matr_tran(4, 4) = (p(1, 1)*p(2, 2)+p(1, 2)*p(2, 1))
        matr_tran(4, 5) = (p(1, 1)*p(2, 3)+p(1, 3)*p(2, 1))
        matr_tran(4, 6) = (p(1, 2)*p(2, 3)+p(1, 3)*p(2, 2))
!
        matr_tran(5, 1) = deux*p(1, 1)*p(3, 1)
        matr_tran(5, 2) = deux*p(1, 2)*p(3, 2)
        matr_tran(5, 3) = deux*p(1, 3)*p(3, 3)
        matr_tran(5, 4) = p(1, 1)*p(3, 2)+p(1, 2)*p(3, 1)
        matr_tran(5, 5) = p(1, 1)*p(3, 3)+p(1, 3)*p(3, 1)
        matr_tran(5, 6) = p(1, 2)*p(3, 3)+p(1, 3)*p(3, 2)
!
        matr_tran(6, 1) = deux*p(2, 1)*p(3, 1)
        matr_tran(6, 2) = deux*p(2, 2)*p(3, 2)
        matr_tran(6, 3) = deux*p(2, 3)*p(3, 3)
        matr_tran(6, 4) = p(2, 1)*p(3, 2)+p(2, 2)*p(3, 1)
        matr_tran(6, 5) = p(2, 1)*p(3, 3)+p(2, 3)*p(3, 1)
        matr_tran(6, 6) = p(2, 2)*p(3, 3)+p(3, 2)*p(2, 3)
!
    end if
!
end subroutine
