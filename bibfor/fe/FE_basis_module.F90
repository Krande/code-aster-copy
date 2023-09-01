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
#include "asterfort/elrfno.h"
#include "asterfort/fe_module.h"
#include "blas/ddot.h"
! --------------------------------------------------------------------------------------------------
!
! fe - generic
!
! Module to generate basis function used for fe methods
!
! --------------------------------------------------------------------------------------------------
!
    type FE_basis
!
        character(len=8) :: typema = " "
        integer :: typeEF = EF_LAGRANGE
        integer :: size

! ----- member function
    contains
        procedure, pass :: initCell => init_cell
        procedure, pass :: initFace => init_face
        procedure, pass :: func => feBSCEval
        procedure, pass :: grad => feBSCGradEv
    end type
!
! --------------------------------------------------------------------------------------------------
    public  :: FE_basis
    private :: feBSCEval, feBSCGradEv, init_cell, init_face
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine init_cell(this, FECell)
!
        implicit none
!
        class(fe_basis), intent(inout) :: this
        type(FE_Cell), intent(in)      :: FECEll
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   Initialization
! --------------------------------------------------------------------------------------------------
!
        this%typema = FECEll%typemas
        this%typeEF = EF_LAGRANGE
!
        if (this%typeEF == EF_LAGRANGE) then
            call elrfno(this%typema, this%size)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine init_face(this, FEFace)
!
        implicit none
!
        class(fe_basis), intent(inout) :: this
        type(FE_Face), intent(in)      :: FEFace
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   Initialization
! --------------------------------------------------------------------------------------------------
!
        this%typema = FEFace%typemas
        this%typeEF = EF_LAGRANGE
!
        if (this%typeEF == EF_LAGRANGE) then
            call elrfno(this%typema, this%size)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function feBSCEval(this, point) result(basisScalEval)
!
        implicit none
!
        class(fe_basis), intent(in)             :: this
        real(kind=8), dimension(3), intent(in)  :: point
        real(kind=8), dimension(MAX_BS)         :: basisScalEval
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   evaluate fe basis scalar
!   In this                 : fe_basis_scalar_cell
!   In point                : point where evaluate
!   Out basisScalEval       : evaluation of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        basisScalEval = 0.d0
        if (this%typeEF == EF_LAGRANGE) then
            call elrfvf(this%typema, point, basisScalEval)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!
!===================================================================================================
!
!===================================================================================================
!
    function feBSCGradEv(this, point) result(BSGradEval)
!
        implicit none
!
        class(fe_basis), intent(in)             :: this
        real(kind=8), dimension(3), intent(in)  :: point
        real(kind=8), dimension(3, MAX_BS)      :: BSGradEval
!
! --------------------------------------------------------------------------------------------------
!   fe - basis functions
!
!   evaluate fe basis scalar
!   In this                 : fe_basis_scalar_cell
!   In point                : point where evaluate
!   Out BSGradEval   : evaluation of the gradient of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        BSGradEval = 0.d0
        if (this%typeEF == EF_LAGRANGE) then
            call elrfdf(this%typema, point, BSGradEval)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
end module
