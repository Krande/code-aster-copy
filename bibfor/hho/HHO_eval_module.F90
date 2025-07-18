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
module HHO_eval_module
!
    use HHO_type
    use HHO_basis_module
    use HHO_quadrature_module
    use HHO_utils_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/elrfvf.h"
#include "asterfort/HHO_size_module.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Evalutaion
!
! Some common utilitaries to evalutaion HHO fonction
!
! --------------------------------------------------------------------------------------------------
    public :: hhoEvalScalCell, hhoEvalScalFace, hhoEvalVecCell, hhoEvalVecFace
    public :: hhoEvalMatCell, hhoEvalSymMatCell, hhoFuncRScalEvalCellQp
    public :: hhoFuncFScalEvalQp, hhoFuncRScalEvalQp, hhoFuncRVecEvalQp, hhoFuncRVecEvalCellQp
!    private  ::
!
contains
!
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoEvalScalCell(hhoBasisCell, order, pt, coeff, size_coeff) result(eval)
!
        implicit none
!
        type(HHO_basis_cell), intent(inout) :: hhoBasisCell
        integer(kind=8), intent(in) :: order
        real(kind=8), dimension(3), intent(in) :: pt
        real(kind=8), dimension(:), intent(in) :: coeff
        integer(kind=8), intent(in) :: size_coeff
        real(kind=8) :: eval
!
! --------------------------------------------------------------------------------------------------
!
!   evaluate a scalar function at a point pt
!   In hhoBasisCell : basis cell
!   In Order        : polynomial order of the function
!   In pt           : point where evaluate
!   In coeff        : polynomial coefficient of the function
!   In size_coeff   : number of coefficient
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: BSCEval
        blas_int :: b_incx, b_incy, b_n
!
        eval = 0.d0
!
! --- Evaluate basis function at pt
        call hhoBasisCell%BSEval(pt, 0, order, BSCEval)
!
        b_n = to_blas_int(size_coeff)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        eval = ddot(b_n, coeff, b_incx, BSCEval, b_incy)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoEvalScalFace(hhoBasisFace, order, pt, coeff, size_coeff) result(eval)
!
        implicit none
!
        type(HHO_basis_face), intent(inout) :: hhoBasisFace
        integer(kind=8), intent(in) :: order
        real(kind=8), dimension(3), intent(in) :: pt
        real(kind=8), dimension(:), intent(in) :: coeff
        integer(kind=8), intent(in) :: size_coeff
        real(kind=8) :: eval
!
! --------------------------------------------------------------------------------------------------
!
!   evaluate a scalar at a point pt
!   In hhoBasisFace : basis Face
!   In Order        : polynomial order of the function
!   In pt           : point where evaluate
!   In coeff        : polynomial coefficient of the function
!   In size_coeff   : number of coefficient
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: BSFEval
        blas_int :: b_incx, b_incy, b_n
!
        eval = 0.d0
!
! --- Evaluate basis function at pt
        call hhoBasisFace%BSEval(pt, 0, order, BSFEval)
!
        b_n = to_blas_int(size_coeff)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        eval = ddot(b_n, coeff, b_incx, BSFEval, b_incy)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoEvalVecCell(hhoBasisCell, order, pt, coeff, size_coeff) result(eval)
!
        implicit none
!
        type(HHO_basis_cell), intent(inout) :: hhoBasisCell
        integer(kind=8), intent(in) :: order
        real(kind=8), dimension(3), intent(in) :: pt
        real(kind=8), dimension(:), intent(in) :: coeff
        integer(kind=8), intent(in) :: size_coeff
        real(kind=8) :: eval(3)
!
! --------------------------------------------------------------------------------------------------
!
!   evaluate a vector at a point pt
!   In hhoBasisCell : basis cell
!   In Order        : polynomial order of the function
!   In pt           : point where evaluate
!   In coeff        : polynomial coefficient of the function
!   In size_coeff   : number of coefficient
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: BSCEval
        integer(kind=8) :: i, size_cmp, deca
        blas_int :: b_incx, b_incy, b_n
!
        eval = 0.d0
        size_cmp = size_coeff/hhoBasisCell%ndim
!
! --- Evaluate basis function at pt
        call hhoBasisCell%BSEval(pt, 0, order, BSCEval)
!
        deca = 0
        do i = 1, hhoBasisCell%ndim
            b_n = to_blas_int(size_cmp)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            eval(i) = ddot(b_n, coeff(deca+1:deca+size_cmp), b_incx, BSCEval, b_incy)
            deca = deca+size_cmp
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoEvalVecFace(hhoBasisFace, order, pt, coeff, size_coeff) result(eval)
!
        implicit none
!
        type(HHO_basis_face), intent(inout) :: hhoBasisFace
        integer(kind=8), intent(in) :: order
        real(kind=8), dimension(3), intent(in) :: pt
        real(kind=8), dimension(:), intent(in) :: coeff
        integer(kind=8), intent(in) :: size_coeff
        real(kind=8) :: eval(3)
!
! --------------------------------------------------------------------------------------------------
!
!   evaluate a vector at a point pt
!   In hhoBasisFace : basis Face
!   In Order        : polynomial order of the function
!   In pt           : point where evaluate
!   In coeff        : polynomial coefficient of the function
!   In size_coeff   : number of coefficient
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: BSFEval
        integer(kind=8) :: i, size_cmp, deca
        blas_int :: b_incx, b_incy, b_n
!
        eval = 0.d0
        size_cmp = size_coeff/(hhoBasisFace%ndim+1)
!
! --- Evaluate basis function at pt
        call hhoBasisFace%BSEval(pt, 0, order, BSFEval)
!
        deca = 0
        do i = 1, (hhoBasisFace%ndim+1)
            b_n = to_blas_int(size_cmp)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            eval(i) = ddot(b_n, coeff(deca+1:deca+size_cmp), b_incx, BSFEval, b_incy)
            deca = deca+size_cmp
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoEvalMatCell(hhoBasisCell, order, pt, coeff, size_coeff) result(eval)
!
        implicit none
!
        type(HHO_basis_cell), intent(inout) :: hhoBasisCell
        integer(kind=8), intent(in) :: order
        real(kind=8), dimension(3), intent(in) :: pt
        real(kind=8), dimension(:), intent(in) :: coeff
        integer(kind=8), intent(in) :: size_coeff
        real(kind=8) :: eval(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!   evaluate a matrix at a point pt
!   In hhoBasisCell : basis cell
!   In Order        : polynomial order of the function
!   In pt           : point where evaluate
!   In coeff        : polynomial coefficient of the function
!   In size_coeff   : number of coefficient
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: BSCEval
        integer(kind=8) :: i, j, size_cmp, deca
        blas_int :: b_incx, b_incy, b_n
!
        eval = 0.d0
        size_cmp = size_coeff/(hhoBasisCell%ndim*hhoBasisCell%ndim)
!
! --- Evaluate basis function at pt
        call hhoBasisCell%BSEval(pt, 0, order, BSCEval)
!
        deca = 0
        do i = 1, hhoBasisCell%ndim
            do j = 1, hhoBasisCell%ndim
                b_n = to_blas_int(size_cmp)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                eval(i, j) = ddot(b_n, coeff(deca+1:deca+size_cmp), b_incx, BSCEval, b_incy)
                deca = deca+size_cmp
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoEvalSymMatCell(hhoBasisCell, order, pt, coeff, size_coeff) result(eval)
!
        implicit none
!
        type(HHO_basis_cell), intent(inout) :: hhoBasisCell
        integer(kind=8), intent(in) :: order
        real(kind=8), dimension(3), intent(in) :: pt
        real(kind=8), dimension(:), intent(in) :: coeff
        integer(kind=8), intent(in) :: size_coeff
        real(kind=8) :: eval(6)
!
! --------------------------------------------------------------------------------------------------
!
!   evaluate a symetrix matrix at a point pt
!   In hhoBasisCell : basis cell
!   In Order        : polynomial order of the function
!   In pt           : point where evaluate
!   In coeff        : polynomial coefficient of the function
!   In size_coeff   : number of coefficient
!   Out mat         : format (XX YY ZZ SQRT(2)*XY SQRT(2)*XZ SQRT(2)*YZ)
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: BSCEval
        real(kind=8) :: mat(3, 3)
        integer(kind=8) :: i, j, size_cmp, deca
        blas_int :: b_incx, b_incy, b_n
!
        if (hhoBasisCell%ndim == 2) then
            size_cmp = size_coeff/3
        else if (hhoBasisCell%ndim == 3) then
            size_cmp = size_coeff/6
        else
            ASSERT(ASTER_FALSE)
        end if
!
! --- Evaluate basis function at pt
        call hhoBasisCell%BSEval(pt, 0, order, BSCEval)
!
        deca = 0
        mat = 0.d0
        do i = 1, hhoBasisCell%ndim
            b_n = to_blas_int(size_cmp)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            mat(i, i) = ddot(b_n, coeff(deca+1:deca+size_cmp), b_incx, BSCEval, b_incy)
            deca = deca+size_cmp
        end do
!
        do i = 1, hhoBasisCell%ndim
            do j = i+1, hhoBasisCell%ndim
                b_n = to_blas_int(size_cmp)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                mat(i, j) = ddot(b_n, coeff(deca+1:deca+size_cmp), b_incx, BSCEval, b_incy)
                mat(j, i) = mat(i, j)
                deca = deca+size_cmp
            end do
        end do
!
        ASSERT(deca == size_coeff)
!
        eval = 0.d0
!
        eval(1) = mat(1, 1)
        eval(2) = mat(2, 2)
        eval(3) = mat(3, 3)
! ---- We don't multiply extra-digonal terms by srqt(2) since we have already multiply
! ---- the basis function by srqt(2)
        eval(4) = mat(1, 2)
        eval(5) = mat(1, 3)
        eval(6) = mat(2, 3)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoFuncFScalEvalQp(hhoQuad, nomfunc, nbpara, nompara, valpara, &
                                  ndim, FuncValuesQp, coeff_mult)
!
        implicit none
!
        type(HHO_Quadrature), intent(in) :: hhoQuad
        character(len=8), intent(in) :: nomfunc
        integer(kind=8), intent(in) :: nbpara
        character(len=8), intent(in) :: nompara(*)
        real(kind=8), intent(inout) :: valpara(*)
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(out) :: FuncValuesQP(*)
        real(kind=8), optional, intent(in) :: coeff_mult
!
!
! --------------------------------------------------------------------------------------------------
!   HHO - Evaluation
!
!   Evaluate an analytical function (*_F) at the quadrature points F(X, Y, ...)
!   In hhoQuad  : Quadrature
!   In nomfunc  : name of the function
!   In nbpara   : number of parameter of the function
!   In nompara  : name of parameters
!   In valpara  : values of parameter
!   In ndim     : spacial dimension of the function (the coordinate-parameters (X,Y,Z) are
!                 always the first parameters) (ndim=0, if the function does not depend on (X,Y,Z))
!   Out FuncValues : values of the function at the quadrature points
!   In coeff_mult  : multply all values by this coefficient (optional)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: npg, ipg, iret
        blas_int :: b_incx, b_n
!
        npg = hhoQuad%nbQuadPoints
! ---- Value of the function at the quadrature point
!
        if (ndim == 0) then
            do ipg = 1, npg
                call fointe('FM', nomfunc, nbpara, nompara, valpara, &
                            FuncValuesQP(ipg), iret)
                ASSERT(iret == 0)
            end do
        else if (ndim <= 3) then
            do ipg = 1, npg
                valpara(1:ndim) = hhoQuad%points(1:ndim, ipg)
                call fointe('FM', nomfunc, nbpara, nompara, valpara, &
                            FuncValuesQP(ipg), iret)
                ASSERT(iret == 0)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (present(coeff_mult)) then
            b_n = to_blas_int(npg)
            b_incx = to_blas_int(1)
            call dscal(b_n, coeff_mult, FuncValuesQP, b_incx)
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoFuncRScalEvalQp(hhoFace, hhoQuad, funcnoEF, FuncValuesQp, coeff_mult)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: funcnoEF(*)
        real(kind=8), intent(out) :: FuncValuesQP(MAX_QP_FACE)
        real(kind=8), optional, intent(in) :: coeff_mult
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Evaluate a function (*_R) at the quadrature points (given at the nodes)
!   In hhoQuad  : Quadrature
!   In funcnoEF : values of the function at the nodes of the EF face
!   Out FuncValues : values of the function at the quadrature points
!   In coeff_mult  : multply all values by this coefficient (optional)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: npg, ipg, ino
        real(kind=8) :: ff(9)
        character(len=8) :: typma
        blas_int :: b_incx, b_n
!
        call cellNameL2S(hhoFace%typema, typma)
        FuncValuesQP = 0.d0
        npg = hhoQuad%nbQuadPoints
        ASSERT(npg <= MAX_QP_FACE)
        ASSERT(hhoQuad%l_point_param)
!
        ff = 0.d0
        do ipg = 1, npg
            call elrfvf(typma, hhoQuad%points_param(1:3, ipg), ff)
            do ino = 1, hhoFace%nbnodes
                FuncValuesQP(ipg) = FuncValuesQP(ipg)+ff(ino)*funcnoEF(ino)
            end do
        end do
!
        if (present(coeff_mult)) then
            b_n = to_blas_int(npg)
            b_incx = to_blas_int(1)
            call dscal(b_n, coeff_mult, FuncValuesQP, b_incx)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoFuncRVecEvalQp(hhoFace, hhoQuad, funcnoEF, FuncValuesQp, coeff_mult)
!
        implicit none
!
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: funcnoEF(*)
        real(kind=8), intent(out) :: FuncValuesQP(3, MAX_QP_FACE)
        real(kind=8), optional, intent(in) :: coeff_mult
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Evaluate a function (*_R) at the quadrature points (given at the nodes)
!   In hhoFace  : Face HHO (!!! The EF face have to be plane)
!   In hhoQuad  : Quadrature
!   In funcnoEF : values of the function at the nodes of the EF face
!   Out FuncValues : values of the function at the quadrature points
!   In coeff_mult  : multply all values by this coefficient (optional)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: npg, ino, idim, celldim, ipg
        real(kind=8) :: ff(9)
        character(len=8) :: typma
        blas_int :: b_incx, b_n
!
        FuncValuesQP = 0.d0
        npg = hhoQuad%nbQuadPoints
        celldim = hhoFace%ndim+1
        ASSERT(npg <= MAX_QP_FACE)
        ASSERT(hhoQuad%l_point_param)
        call cellNameL2S(hhoFace%typema, typma)
!
        do ipg = 1, npg
            call elrfvf(typma, hhoQuad%points_param(1:3, ipg), ff)
            do idim = 1, celldim
                do ino = 1, hhoFace%nbnodes
                    FuncValuesQP(idim, ipg) = FuncValuesQP(idim, ipg)+ff(ino)*funcnoEF(celldim*(&
                                              &ino-1)+idim)
                end do
            end do
        end do
!
        if (present(coeff_mult)) then
            b_n = to_blas_int(3*npg)
            b_incx = to_blas_int(1)
            call dscal(b_n, coeff_mult, FuncValuesQP, b_incx)
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoFuncRVecEvalCellQp(hhoCell, hhoQuad, funcnoEF, FuncValuesQp, coeff_mult)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: funcnoEF(*)
        real(kind=8), intent(out) :: FuncValuesQP(3, MAX_QP_CELL)
        real(kind=8), optional, intent(in) :: coeff_mult
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Evaluate a function (*_R) at the quadrature points (given at the nodes)
!   In hhoCell  : Cell HHO
!   In hhoQuad  : Quadrature
!   In funcnoEF : values of the function at the nodes of the EF cell
!   Out FuncValues : values of the function at the quadrature points
!   In coeff_mult  : multply all values by this coefficient (optional)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: npg, ino, idim, ipg
        real(kind=8) :: ff(27)
        character(len=8) :: typma
        blas_int :: b_incx, b_n
!
        FuncValuesQP = 0.d0
        npg = hhoQuad%nbQuadPoints
        ASSERT(npg <= MAX_QP_CELL)
        ASSERT(hhoQuad%l_point_param)
        call cellNameL2S(hhoCell%typema, typma)
!
        do ipg = 1, npg
            call elrfvf(typma, hhoQuad%points_param(1:3, ipg), ff)
            do idim = 1, hhoCell%ndim
                do ino = 1, hhoCell%nbnodes
                    FuncValuesQP(idim, ipg) = FuncValuesQP(idim, ipg)+ff(ino)*funcnoEF(hhoCell%n&
                                              &dim*(ino-1)+idim)
                end do
            end do
        end do
!
        if (present(coeff_mult)) then
            b_n = to_blas_int(3*npg)
            b_incx = to_blas_int(1)
            call dscal(b_n, coeff_mult, FuncValuesQP, b_incx)
        end if
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoFuncRScalEvalCellQp(hhoCell, hhoQuad, funcnoEF, FuncValuesQp, coeff_mult)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8), intent(in) :: funcnoEF(*)
        real(kind=8), intent(out) :: FuncValuesQP(MAX_QP_CELL)
        real(kind=8), optional, intent(in) :: coeff_mult
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Evaluate a function (*_R) at the quadrature points (given at the nodes)
!   In hhoCell  : Cell HHO
!   In hhoQuad  : Quadrature
!   In funcnoEF : values of the function at the nodes of the EF cell
!   Out FuncValues : values of the function at the quadrature points
!   In coeff_mult  : multply all values by this coefficient (optional)
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: npg, ino, ipg
        real(kind=8) :: ff(27)
        character(len=8) :: typma
        blas_int :: b_incx, b_n
!
        FuncValuesQP = 0.d0
        npg = hhoQuad%nbQuadPoints
        ASSERT(npg <= MAX_QP_CELL)
        ASSERT(hhoQuad%l_point_param)
        call cellNameL2S(hhoCell%typema, typma)
!
        ff = 0.d0
        do ipg = 1, npg
            call elrfvf(typma, hhoQuad%points_param(1:3, ipg), ff)
            do ino = 1, hhoCell%nbnodes
                FuncValuesQP(ipg) = FuncValuesQP(ipg)+ff(ino)*funcnoEF(ino)
            end do
        end do
!
        if (present(coeff_mult)) then
            b_n = to_blas_int(npg)
            b_incx = to_blas_int(1)
            call dscal(b_n, coeff_mult, FuncValuesQP, b_incx)
        end if
!
    end subroutine
!
end module
