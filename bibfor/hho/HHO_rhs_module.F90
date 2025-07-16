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
! person_in_charge: mickael.abbas at edf.fr
!
module HHO_rhs_module
!
    use HHO_basis_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use FE_algebra_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/HHO_size_module.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Module to compute the rhs member (f,v)
!
! --------------------------------------------------------------------------------------------------
!
    public :: hhoMakeRhsFaceScal, hhoMakeRhsFaceVec
    public :: hhoMakeRhsCellScal, hhoMakeRhsCellVec
!    private  ::
!
contains
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsFaceScal(hhoFace, hhoQuad, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: ValuesQP(MAX_QP_FACE)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_FACE_SCAL)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vF)_F
!   In hhoFace      : the current HHO Face
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of scalar function f at the quadrature points
!   In degree       : degree of the basis
!   Out rhs         : (f, vF)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: ipg, size
        type(HHO_basis_face) :: hhoBasisFace
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: BSFEval
        real(kind=8) :: coeff
!
        rhs = 0.d0
!
! -- init face basis
        call hhoBasisFace%initialize(hhoFace)
        size = hhoBasisFace%BSSize(0, degree)
!
! -- Loop on quadrature point
        do ipg = 1, hhoQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), 0, degree, BSFEval)
!
! ---- rhs = rhs + weight * BSFEval
            coeff = hhoQuad%weights(ipg)*ValuesQP(ipg)
            call daxpy_1(size, coeff, BSFEval, rhs)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsFaceVec(hhoFace, hhoQuad, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: ValuesQP(3, MAX_QP_FACE)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_FACE_VEC)
!
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vF)_F
!   In hhoFace      : the current HHO Face
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of vectorial function f at the quadrature points
!   In degree  : degree of the basis
!   Out rhs  : term (f, vF)_F (rhs member)
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_face) :: hhoBasisFace
        integer(kind=8) :: size, idir, begin, i
        real(kind=8) :: Values(MAX_QP_FACE), rhs_dir(MSIZE_FACE_SCAL)
        blas_int :: b_incx, b_incy, b_n
!
! -- init face basis
        call hhoBasisFace%initialize(hhoFace)
        size = hhoBasisFace%BSSize(0, degree)
!
        rhs = 0.d0
        Values = 0.d0
!
        begin = 1
        do idir = 1, hhoFace%ndim+1
            do i = 1, hhoQuad%nbQuadPoints
                Values(i) = ValuesQP(idir, i)
            end do
            call hhoMakeRhsFaceScal(hhoFace, hhoQuad, Values, degree, rhs_dir)
            b_n = to_blas_int(size)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, rhs_dir, b_incx, rhs(begin), b_incy)
            begin = begin+size
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsCellScal(hhoCell, hhoQuad, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: ValuesQP(MAX_QP_CELL)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_CELL_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vT)_T
!   In hhoCell      : the current HHO Cell
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of scalar function f at the quadrature points
!   Out rhs         : (f, vT)_T term
!   In  degree      : degree of the basis
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer(kind=8) :: ipg, size
        type(HHO_basis_cell) :: hhoBasisCell
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: BSCEval
        real(kind=8) :: coeff
!
! -- init face basis
        call hhoBasisCell%initialize(hhoCell)
        size = hhoBasisCell%BSSize(0, degree)
!
        rhs = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, hhoQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, degree, BSCEval)
!
! ---- rhs = rhs + weight * BSCEval
            coeff = hhoQuad%weights(ipg)*ValuesQP(ipg)
            call daxpy_1(size, coeff, BSCEval, rhs)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMakeRhsCellVec(hhoCell, hhoQuad, ValuesQP, degree, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: ValuesQP(3, MAX_QP_CELL)
        integer(kind=8), intent(in) :: degree
        real(kind=8), intent(out) :: rhs(MSIZE_CELL_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the term (f, vT)_T
!   In hhoCell      : the current HHO Cell
!   In hhoQuad      : Quadrature for the face
!   In ValuesQP : Values of vectorial function f at the quadrature points
!   Out rhs         : term (f, vT)_T (rhs member)
!   In degree  : degree of the basis
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        integer(kind=8) :: size, idir, begin, i
        real(kind=8) :: Values(MAX_QP_CELL), rhs_dir(MSIZE_CELL_SCAL)
        blas_int :: b_incx, b_incy, b_n
!
! -- init face basis
        call hhoBasisCell%initialize(hhoCell)
        size = hhoBasisCell%BSSize(0, degree)
!
        rhs = 0.d0
        Values = 0.d0
!
        begin = 1
        do idir = 1, hhoCell%ndim
            do i = 1, hhoQuad%nbQuadPoints
                Values(i) = ValuesQP(idir, i)
            end do
            call hhoMakeRhsCellScal(hhoCell, hhoQuad, Values, degree, rhs_dir)
            b_n = to_blas_int(size)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, rhs_dir, b_incx, rhs(begin), b_incy)
            begin = begin+size
        end do
!
    end subroutine
!
end module
