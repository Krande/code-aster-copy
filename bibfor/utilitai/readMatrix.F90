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
subroutine readMatrix(name, nrows, ncols, l_sym, mat)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
!
    character(len=*), intent(in) :: name
    integer(kind=8), intent(in) :: nrows, ncols
    aster_logical, intent(in) :: l_sym
    real(kind=8), dimension(:, :), intent(inout) :: mat
!
! --------------------------------------------------------------------------------------------------
!
! IO routine
!
! Read a matrix in memory (zr)
!
! --------------------------------------------------------------------------------------------------
!
!   In name        : name of the matrix read with jevech
!   In nrows/ncols : size of the matrix
!   In l_sym       : the matrix is symmetric ?
!   In mat         : matrix to write
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jv_matr_out, j, ij, i
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call jevech(name, 'L', jv_matr_out)
!
    if (l_sym) then
        ASSERT(ncols == nrows)
        do j = 1, ncols
            ij = (j-1)*j/2
            b_n = to_blas_int(j)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jv_matr_out+ij), b_incx, mat(:, j), b_incy)
            do i = 1, j-1
                mat(j, i) = mat(i, j)
            end do
        end do
    else
        do i = 1, nrows
            ij = (i-1)*ncols
            b_n = to_blas_int(ncols)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jv_matr_out+ij), b_incx, mat(i, :), b_incy)
        end do
    end if
!
end subroutine
