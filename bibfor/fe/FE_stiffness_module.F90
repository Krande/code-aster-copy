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
!
    implicit none
!
    private
#include "asterf_types.h"
#include "FE_module.h"
#include "blas/dgemv.h"
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
    public :: FEStiffVecScal, FEStiffMatScal
!    private  ::
!
contains
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
! ----- Local variables
        integer :: ipg
        real(kind=8), dimension(3, MAX_BS) :: BSEval
!
        vec = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%grad(FEQuad%points_param(1:3, ipg))
!
            call dgemv('T', FEQuad%ndim, FEBasis%size, FEQuad%weights(ipg), BSEval, 3, &
                       ValuesQP(1:3, ipg), 1, 1.d0, vec, 1)
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
!   HHO
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
        integer :: ipg, j
        real(kind=8), dimension(3, MAX_BS) :: BSEval
        real(kind=8) :: Kgradj(3)
!
        mat = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            BSEval = FEBasis%grad(FEQuad%points_param(1:3, ipg))
!
            do j = 1, FEBasis%size
                Kgradj = matmul(ValuesQP(1:3, 1:3, ipg), BSEval(1:3, j))
                call dgemv('T', FEQuad%ndim, FEBasis%size, FEQuad%weights(ipg), BSEval, 3, &
                           Kgradj, 1, 1.d0, mat(:, j), 1)
            end do
!
        end do
!
    end subroutine
!
end module
