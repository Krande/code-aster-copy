! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
module HHO_massmat_module
!
    use HHO_basis_module
    use HHO_quadrature_module
    use HHO_type
    use HHO_utils_module
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/HHO_size_module.h"
#include "blas/dsyr.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Module to compute mass matrix
!
! --------------------------------------------------------------------------------------------------
    type HHO_massmat_cell
        integer :: max_nrows = MSIZE_CELL_SCAL, max_ncols = MSIZE_CELL_SCAL
        integer :: nrows = 0, ncols = 0
        aster_logical :: isIdentity = ASTER_FALSE
        real(kind=8) :: m(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)

! ----- member function
    contains
        procedure, pass :: compute => hhoMassMatCellScal

    end type
!
    type HHO_massmat_face
        integer :: max_nrows = MSIZE_FACE_SCAL, max_ncols = MSIZE_FACE_SCAL
        integer :: nrows = 0, ncols = 0
        aster_logical :: isIdentity = ASTER_FALSE
        real(kind=8) :: m(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL)

! ----- member function
    contains
        procedure, pass :: compute => hhoMassMatFaceScal

    end type
!
    public :: HHO_massmat_cell, HHO_massmat_face
    private  :: hhoMassMatCellScal, hhoMassMatFaceScal
!    private  ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMassMatCellScal(this, hhoCell, min_order, max_order)
!
        implicit none
!
        class(HHO_massmat_cell), intent(inout) :: this
        type(HHO_Cell), intent(in) :: hhoCell
        integer, intent(in)        :: min_order
        integer, intent(in)        :: max_order
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the scalar mass matrix of a Cell form order "min_order" to order "max_order"
!   In hhoCell     : the current HHO Cell
!   In min_order   : minimum order to evaluate
!   In max_order   : maximum order to evaluate
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_quadrature)  :: hhoQuad
        real(kind=8), dimension(MSIZE_CELL_SCAL):: basisScalEval
        integer :: dimMat, ipg, i
        aster_logical :: dbg
! --------------------------------------------------------------------------------------------------
!
! ----- init basis
        call hhoBasisCell%initialize(hhoCell)
! ----- dimension of massMat
        dimMat = hhoBasisCell%BSSize(min_order, max_order)
        this%m = 0.d0
        this%nrows = dimMat
        this%ncols = dimMat
!
#ifdef ASTER_DEBUG_FC
        dbg = ASTER_TRUE
#else
        dbg = ASTER_FALSE
#endif
!
        if (hhoBasisCell%isOrthonormal() .and. .not. dbg .and. &
            (hhoCell%typema .eq. 'TETRA4' .or. hhoCell%ndim .eq. 2)) then
            do i = 1, dimMat
                this%m(i, i) = 1.d0
            end do
            this%isIdentity = ASTER_TRUE
        else
!
! ----- get quadrature
            call hhoQuad%GetQuadCell(hhoCell, 2*max_order)
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval bais function at the quadrature point
                call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), min_order, max_order,&
                                  & basisScalEval)
! --------  Eval massMat
                call dsyr('U', dimMat, hhoQuad%weights(ipg), basisScalEval, &
                          1, this%m, this%max_nrows)
            end do
!
! ----- Copy the lower part
!
            call hhoCopySymPartMat('U', this%m(1:dimMat, 1:dimMat))
!
            if (hhoBasisCell%isOrthonormal()) then
                this%isIdentity = hhoIsIdentityMat(this%m, dimMat)
#ifdef ASTER_DEBUG_FC
                ASSERT(this%isIdentity)
#endif
            end if

        end if
        ! call hhoPrintMat(this%m(1:dimMat, 1:dimMat))
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMassMatFaceScal(this, hhoFace, min_order, max_order)
!
        implicit none
!
        class(HHO_massmat_face), intent(inout) :: this
        type(HHO_Face), intent(in) :: hhoFace
        integer, intent(in)        :: min_order
        integer, intent(in)        :: max_order
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the scalar mass matrix of a Face form order "min_order" to order "max_order"
!   In hhoFace     : the current HHO Face
!   In min_order   : minimum order to evaluate
!   In max_order   : maximum order to evaluate
!
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_face) :: hhoBasisFace
        type(HHO_quadrature)  :: hhoQuad
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: basisScalEval
        integer :: dimMat, ipg, i
        aster_logical :: dbg
! --------------------------------------------------------------------------------------------------
!
! ----- init basis
        call hhoBasisFace%initialize(hhoFace)
! ----- dimension of massMat
        dimMat = hhoBasisFace%BSSize(min_order, max_order)
!
        this%m = 0.d0
        this%nrows = dimMat
        this%ncols = dimMat
!
#ifdef ASTER_DEBUG_FC
        dbg = ASTER_TRUE
#else
        dbg = ASTER_FALSE
#endif
!
        if (hhoBasisFace%isOrthonormal() .and. .not. dbg) then
            do i = 1, dimMat
                this%m(i, i) = 1.d0
            end do
            this%isIdentity = ASTER_TRUE
        else
!
! ----- get quadrature
            call hhoQuad%GetQuadFace(hhoFace, 2*max_order)
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval bais function at the quadrature point
                call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), min_order, &
                                        & max_order, basisScalEval)
! --------  Eval massMat
                call dsyr('U', dimMat, hhoQuad%weights(ipg), basisScalEval, &
                          1, this%m, this%max_nrows)
            end do
!
! ----- Copy the lower part
!
            call hhoCopySymPartMat('U', this%m(1:dimMat, 1:dimMat))
!
            if (hhoBasisFace%isOrthonormal()) then
                this%isIdentity = hhoIsIdentityMat(this%m, dimMat)
#ifdef ASTER_DEBUG_FC
                ASSERT(this%isIdentity)
#endif
            end if
        end if
!
        ! call hhoPrintMat(this%m)
!
    end subroutine
!
end module
