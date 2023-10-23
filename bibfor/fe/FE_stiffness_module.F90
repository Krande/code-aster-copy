! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
module FE_stiffness_module
!
    use FE_basis_module
    use FE_quadrature_module
    use HHO_utils_module, only: hhoCopySymPartMat
!
    implicit none
!
    private
#include "asterf_types.h"
#include "FE_module.h"
#include "blas/dgemv.h"
#include "asterfort/assert.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE - Finite Element
!
! Module to compute stiffness terms
!
! --------------------------------------------------------------------------------------------------
!
    public :: FEStiffVecScal, FEStiffMatScal, FEStiffVecScalAdd, FEStiffMatScalAdd
!    private  ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffVecScalAdd(FEBasis, BGSEval, weight, ValuesQP, vec)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(3, MAX_BS)  :: BGSEval
        real(kind=8), intent(in)            :: weight
        real(kind=8), intent(inout)         :: vec(MAX_BS)
        real(kind=8), intent(in)            :: ValuesQP(3)
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the rigidity vector
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
!
!
        call dgemv('T', FEBasis%ndim, FEBasis%size, weight, BGSEval, 3, &
                   ValuesQP, 1, 1.d0, vec, 1)
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffVecScal(FEQuad, FEBasis, ValuesQP, vec)
!
        implicit none
!
        type(FE_Quadrature), intent(in)     :: FEQuad
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(out)           :: vec(MAX_BS)
        real(kind=8), intent(in)            :: ValuesQP(3, MAX_QP)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity vector
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer :: ipg
        real(kind=8), dimension(3, MAX_BS) :: BSEval
!
        vec = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%grad(FEQuad%points_param(1:3, ipg), FEQuad%jacob(1:3, 1:3, ipg))
!
            call FEStiffVecScalAdd(FEBasis, BSEval, FEQuad%weights(ipg), ValuesQP(1:3, ipg), vec)
        end do
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffMatScalAdd(FEBasis, BGSEval, weight, ValueQP, mat)
!
        implicit none
!
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(in), dimension(3, MAX_BS) :: BGSEval
        real(kind=8), intent(in)            :: weight
        real(kind=8), intent(inout)           :: mat(MAX_BS, MAX_BS)
        real(kind=8), intent(in)            :: ValueQP(3, 3)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity matrix
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer :: j
        real(kind=8) :: Kgradj(3)
!
        do j = 1, FEBasis%size
            if (FEBasis%ndim == 3) then
                Kgradj(1) = ValueQP(1, 1)*BGSEval(1, j)+ValueQP(1, 2)*BGSEval(2, j)+ &
                            ValueQP(1, 3)*BGSEval(3, j)
                Kgradj(2) = ValueQP(2, 1)*BGSEval(1, j)+ValueQP(2, 2)*BGSEval(2, j)+ &
                            ValueQP(2, 3)*BGSEval(3, j)
                Kgradj(3) = ValueQP(3, 1)*BGSEval(1, j)+ValueQP(3, 2)*BGSEval(2, j)+ &
                            ValueQP(3, 3)*BGSEval(3, j)
            elseif (FEBasis%ndim == 2) then
                Kgradj(1) = ValueQP(1, 1)*BGSEval(1, j)+ValueQP(1, 2)*BGSEval(2, j)
                Kgradj(2) = ValueQP(2, 1)*BGSEval(1, j)+ValueQP(2, 2)*BGSEval(2, j)
            else
                Kgradj(1) = ValueQP(1, 1)*BGSEval(1, j)
            end if
            call dgemv('T', FEBasis%ndim, j, weight, BGSEval, 3, Kgradj, 1, 1.d0, mat(:, j), 1)
        end do
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEStiffMatScal(FEQuad, FEBasis, ValuesQP, mat)
!
        implicit none
!
        type(FE_Quadrature), intent(in)     :: FEQuad
        type(FE_Basis), intent(in)          :: FEBasis
        real(kind=8), intent(out)           :: mat(MAX_BS, MAX_BS)
        real(kind=8), intent(in)            :: ValuesQP(3, 3, MAX_QP)
! --------------------------------------------------------------------------------------------------
!
!
!   Compute the rigidity matrix
!   In hhoQuad      : Quadrature
!   In hhoBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, grad v)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer :: ipg
        real(kind=8), dimension(3, MAX_BS) :: BSEval
!
        mat = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%grad(FEQuad%points_param(1:3, ipg), FEQuad%jacob(1:3, 1:3, ipg))

            call FEStiffMatScalAdd(FEBasis, BSEval, FEQuad%weights(ipg), ValuesQP(1:3, 1:3, ipg), &
                                   mat)
!
        end do
!
! ----- Copy the lower part
!
        call hhoCopySymPartMat('U', mat(1:FEBasis%size, 1:FEBasis%size))
!
    end subroutine
!
end module
