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
subroutine readVector(name, nb_value, vect, offset)
!
    implicit none
!
#include "asterfort/jevech.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
!
    character(len=*), intent(in) :: name
    integer(kind=8), intent(in) :: nb_value
    real(kind=8), intent(out) :: vect(*)
    integer(kind=8), intent(in), optional :: offset
!
! --------------------------------------------------------------------------------------------------
!
! IO routine
!
! Read a vector in memory (zr)
!
! --------------------------------------------------------------------------------------------------
!
!   In name        : name of the vector read with jevech
!   In nb_value    : size of the vector
!   Out vect       : vector to read
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jv_vect_in, offset_
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    offset_ = 0
    if (present(offset)) offset_ = offset
!
    call jevech(name, 'L', jv_vect_in)
    b_n = to_blas_int(nb_value)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jv_vect_in+offset_), b_incx, vect, b_incy)
!
end subroutine
