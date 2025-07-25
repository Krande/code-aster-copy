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
module bloc_fe_module
!
! Utilitaires pour faciliter la manipulation d'un element fini par bloc (ex: EF mixte)
! Encapsule l'appel aux blas
!
    implicit none
!
#include "asterfort/assert.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
!
    public :: prod_bd
    public :: prod_sb
    public :: prod_bkb
    public :: add_fint
    public :: add_matr
!
contains
!
    function prod_bd(b, d) result(bd)
        implicit none
        real(kind=8) :: b(:, :, :), d(:, :)
        real(kind=8) :: bd(size(b, 1))
! -------------------------------------------------------------------------
        integer(kind=8) :: neps, ndim, nno
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
! -------------------------------------------------------------------------
        neps = size(b, 1)
        ndim = size(b, 2)
        nno = size(b, 3)
        ASSERT(size(d, 1) .eq. ndim)
        ASSERT(size(d, 2) .eq. nno)
        b_lda = to_blas_int(neps)
        b_m = to_blas_int(neps)
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('n', b_m, b_n, 1.d0, b, &
                   b_lda, d, b_incx, 0.d0, bd, &
                   b_incy)
    end function prod_bd
!
    function prod_sb(s, b) result(sb)
        implicit none
        real(kind=8) :: s(:), b(:, :, :)
        real(kind=8) :: sb(size(b, 2), size(b, 3))
! -------------------------------------------------------------------------
        integer(kind=8) :: ndim, nno, neu
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
! -------------------------------------------------------------------------
        ndim = size(b, 2)
        nno = size(b, 3)
        neu = size(s)
        ASSERT(neu .eq. size(b, 1))
        b_lda = to_blas_int(neu)
        b_m = to_blas_int(neu)
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('t', b_m, b_n, 1.d0, b, &
                   b_lda, s, b_incx, 0.d0, sb, &
                   b_incy)
    end function prod_sb
!
    function prod_bkb(bin, kinjm, bjm) result(bkb)
        implicit none
        real(kind=8), intent(in) :: bin(:, :, :), bjm(:, :, :), kinjm(:, :)
        real(kind=8) :: bkb(size(bin, 2), size(bin, 3), size(bjm, 2), size(bjm, 3))
! -------------------------------------------------------------------------
        integer(kind=8) :: ndim1, ndim2, nno1, nno2, neps1, neps2
        real(kind=8), allocatable :: kbjm(:, :, :)
        blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m, b_n
! -------------------------------------------------------------------------
        neps1 = size(bin, 1)
        ndim1 = size(bin, 2)
        nno1 = size(bin, 3)
        neps2 = size(bjm, 1)
        ndim2 = size(bjm, 2)
        nno2 = size(bjm, 3)
        ASSERT(size(kinjm, 1) .eq. neps1)
        ASSERT(size(kinjm, 2) .eq. neps2)
!
        allocate (kbjm(neps1, ndim2, nno2))
!
        b_ldc = to_blas_int(neps1)
        b_ldb = to_blas_int(neps2)
        b_lda = to_blas_int(neps1)
        b_m = to_blas_int(neps1)
        b_n = to_blas_int(ndim2*nno2)
        b_k = to_blas_int(neps2)
        call dgemm('n', 'n', b_m, b_n, b_k, &
                   1.d0, kinjm, b_lda, bjm, b_ldb, &
                   0.d0, kbjm, b_ldc)
        b_ldc = to_blas_int(ndim1*nno1)
        b_ldb = to_blas_int(neps1)
        b_lda = to_blas_int(neps1)
        b_m = to_blas_int(ndim1*nno1)
        b_n = to_blas_int(ndim2*nno2)
        b_k = to_blas_int(neps1)
        call dgemm('t', 'n', b_m, b_n, b_k, &
                   1.d0, bin, b_lda, kbjm, b_ldb, &
                   0.d0, bkb, b_ldc)
!
        deallocate (kbjm)
    end function prod_bkb
!
    subroutine add_fint(fint, xin, fin)
! Contribution d'un bloc de forces au vecteur des forces interieures
        implicit none
        integer(kind=8), intent(in) :: xin(:, :)
        real(kind=8), intent(in) :: fin(:, :)
        real(kind=8), intent(inout) :: fint(:)
! -------------------------------------------------------------------------
        integer(kind=8) :: ndim, nno, i, n
! -------------------------------------------------------------------------
        ndim = size(fin, 1)
        nno = size(fin, 2)
        ASSERT(ndim .eq. size(xin, 1))
        ASSERT(nno .eq. size(xin, 2))
        forall (i=1:ndim, n=1:nno) fint(xin(i, n)) = fint(xin(i, n))+fin(i, n)
    end subroutine add_fint
!
    subroutine add_matr(matr, xin, xjm, kinjm)
! Contribution d'un bloc de matrice a la matrice tangente (stockage j,m,i,n)
        implicit none
        integer(kind=8), intent(in) :: xin(:, :), xjm(:, :)
        real(kind=8), intent(in) :: kinjm(:, :, :, :)
        real(kind=8), intent(inout) :: matr(:, :)
! -------------------------------------------------------------------------
        integer(kind=8) :: ndimin, nnoin, ndimjm, nnojm, i, n, j, m
! -------------------------------------------------------------------------
        ndimin = size(kinjm, 1)
        nnoin = size(kinjm, 2)
        ndimjm = size(kinjm, 3)
        nnojm = size(kinjm, 4)
        ASSERT(ndimin .eq. size(xin, 1))
        ASSERT(nnoin .eq. size(xin, 2))
        ASSERT(ndimjm .eq. size(xjm, 1))
        ASSERT(nnojm .eq. size(xjm, 2))
        forall (i=1:ndimin, n=1:nnoin, j=1:ndimjm, m=1:nnojm) matr(xjm(j, m), xin(i, n)) = &
            matr(xjm(j, m), xin(i, n))+kinjm(i, n, j, m)
    end subroutine add_matr
!
end module bloc_fe_module
