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
module FE_eval_module
!
    use FE_basis_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "FE_module.h"
#include "blas/ddot.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE - Finite Element
!
! Module to eval function
!
! --------------------------------------------------------------------------------------------------
!
    public :: FEEvalFuncScal, FEEvalGradVec
!    private  ::
!
contains
!
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalFuncScal( FEBasis, val_nodes, point) result(func)
!
        implicit none
!
        type(FE_Basis), intent(in)         :: FEBasis
        real(kind=8), intent(in)           :: val_nodes(*)
        real(kind=8), intent(in)           :: point(3)
        real(kind=8)                       :: func
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate scalar values from value at nodes
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        real(kind=8) :: funcEF(MAX_BS)
!
        funcEF = FEBasis%func(point)
        func =  ddot(FEBasis%size, val_nodes, 1, funcEF, 1)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function FEEvalGradVec( FEBasis, val_nodes, point) result(grad)
!
        implicit none
!
        type(FE_Basis), intent(in)         :: FEBasis
        real(kind=8), intent(in)           :: val_nodes(*)
        real(kind=8), intent(in)           :: point(3)
        real(kind=8)                       :: grad(3)
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Evaluate scalar values from value at nodes
!   In FEBasis     : tBasis function
!   In ValuesQP     : Values of scalar function f at the quadrature points
!   Out rhs         : (f, v)_F term
!
! --------------------------------------------------------------------------------------------------
!
! ----- Local variables
        integer :: i
        real(kind=8) :: gradEF(3, MAX_BS)
!
        grad = 0.d0
        gradEF = FEBasis%grad(point)
        do i = 1, FEBasis%size
            grad = grad+val_nodes(i)*gradEF(1:3, i)
        end do
!
    end function
!
end module
