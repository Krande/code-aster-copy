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
! person_in_charge: mickael.abbas at edf.fr
!
module FE_basis_module
!
    use fe_topo_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrfvf.h"
#include "asterfort/elrfdf.h"
#include "blas/ddot.h"
! --------------------------------------------------------------------------------------------------
!
! fe - generic
!
! Module to generate basis function used for fe methods
!
! --------------------------------------------------------------------------------------------------
!
    type FE_basis_cell

! ----- member function
    contains
        procedure, pass :: func => feBSCEval
        procedure, pass :: grad => feBSCGradEv
    end type
!
    type FE_basis_face

! ----- member function
    contains
        procedure, pass :: func => feBSFEval
    end type
!
! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
    public  :: FE_basis_cell, FE_basis_face
    private :: feBSCEval, feBSCGradEv, feBSFEval
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine feBSFEval(this, FEFace, point, basisScalEval)
!
        implicit none
!
        type(fe_face), intent(in)                              :: feface
        class(fe_basis_face), intent(inout)                    :: this
        real(kind=8), dimension(3), intent(in)                  :: point
        real(kind=8), dimension(9), intent(out)   :: basisScalEval
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   evaluate fe basis scalar for a cell
!   In feCell              : the current fe cell
!   In this                 : fe_basis_scalar_cell
!   In point                : point where evaluate
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out basisScalEval       : evaluation of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        call elrfvf(feface%typemas, point, basisScalEval)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine feBSCEval(this, feCell, point, basisScalEval)
!
        implicit none
!
        type(fe_Cell), intent(in)                              :: feCell
        class(fe_basis_cell), intent(inout)                    :: this
        real(kind=8), dimension(3), intent(in)                  :: point
        real(kind=8), dimension(27), intent(out)   :: basisScalEval
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   evaluate fe basis scalar for a cell
!   In feCell              : the current fe cell
!   In this                 : fe_basis_scalar_cell
!   In point                : point where evaluate
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out basisScalEval       : evaluation of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        call elrfvf(feCell%typemas, point, basisScalEval)
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine feBSCGradEv(this, feCell, point, BSGradEval)
!
        implicit none
!
        type(fe_Cell), intent(in)                              :: feCell
        class(fe_basis_cell), intent(inout)                    :: this
        real(kind=8), dimension(3), intent(in)                  :: point
        real(kind=8), dimension(3, 27), intent(out):: BSGradEval
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   evaluate fe basis scalar for a 3D cell
!   In feCell              : the current fe cell
!   In this                 : fe_basis_scalar_cell
!   In point                : point where evaluate
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out BSGradEval   : evaluation of the gradient of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        call elrfdf(feCell%typemas, point, BSGradEval)

!
    end subroutine
!
end module
