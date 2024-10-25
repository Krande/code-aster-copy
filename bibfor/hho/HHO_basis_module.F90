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
module HHO_basis_module
!
    use HHO_geometry_module
    use HHO_measure_module
    use HHO_monogen_module
    use HHO_quadrature_module
    use HHO_type
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/binomial.h"
#include "asterfort/HHO_basis_module.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/readVector.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/jevech.h"
! --------------------------------------------------------------------------------------------------
!
! HHO - generic
!
! Module to generate basis function used for HHO methods
!
! --------------------------------------------------------------------------------------------------
!
    type HHO_basis_cell
        type(HHO_monomials) :: hhoMono
        integer :: type = BASIS_ORTHO
        integer :: ndim = 0
        real(kind=8) :: center(3) = 0.d0
        real(kind=8) :: scaling_factor(3) = 0.d0
        real(kind=8) :: rotmat(3, 3) = 0.d0
        real(kind=8) :: coeff_mono(MAX_CELL_COEF) = 0.d0
        integer      :: coeff_shift(MSIZE_CELL_SCAL+1) = 0

! ----- member function
    contains
        procedure, pass :: initialize => hhoBasisCellInit
        procedure, pass :: BSSize => hhoBSCellSize
        procedure, pass :: BVSize => hhoBVCellSize
        procedure, pass :: BMSize => hhoBMCellSize
        procedure, pass :: BSRange => hhoBSCellRange
        procedure, pass :: BVRange => hhoBVCellRange
        procedure, pass :: BMRange => hhoBMCellRange
        procedure, pass :: BSEval => hhoBSCellEval
        procedure, pass :: BSEvalGrad => hhoBSCellGradEv
        procedure, pass :: BVEvalSymGrad => hhoBVCellSymGdEv
        procedure, pass :: isOrthonormal => hhoBasisCellType
        procedure, pass, private :: map_pt => map_pt_cell
    end type
!
    type HHO_basis_face
        type(HHO_monomials) :: hhoMono
        integer :: type = BASIS_ORTHO
        integer :: ndim = 0
        real(kind=8) :: center(3) = 0.d0
        real(kind=8) :: scaling_factor(2) = 0.d0
        real(kind=8) :: rotmat(2, 3) = 0.d0
        real(kind=8) :: coeff_mono(MAX_FACE_COEF) = 0.d0
        integer      :: coeff_shift(MSIZE_FACE_SCAL+1) = 0
! ----- member function
    contains
        procedure, pass :: initialize => hhoBasisFaceInit
        procedure, pass :: BSSize => hhoBSFaceSize
        procedure, pass :: BVSize => hhoBVFaceSize
        procedure, pass :: BSRange => hhoBSFaceRange
        procedure, pass :: BVRange => hhoBVFaceRange
        procedure, pass :: BSEval => hhoBSFaceEval
        procedure, pass :: isOrthonormal => hhoBasisFaceType
        procedure, pass, private :: map_pt => map_pt_face
    end type
! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
    public  :: HHO_basis_cell, HHO_basis_face
    private :: hhoBasisCellInit, hhoBasisFaceInit, hhoBSCellSize, hhoBSFaceSize
    private :: hhoBVCellSize, hhoBVFaceSize, hhoBMCellSize
    private :: hhoBSCellRange, hhoBVCellRange, hhoBMCellRange, hhoBSFaceRange, hhoBVFaceRange
    private :: hhoBSCellEval, hhoBSFaceEval, hhoBSCellGradEv, hhoBVCellSymGdEv, check_order
    private :: map_pt_cell, map_pt_face, orthonormalization
    private :: hhoBasisCellType, hhoBasisFaceType
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine check_order(min_order, max_order, maxAutorized)
!
        implicit none
!
        integer, intent(in) :: min_order
        integer, intent(in) :: max_order
        integer, intent(in) :: maxAutorized
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   control the order of integration of [min_order, max_order] has to be included in
!   [0, maxAutorized]
!   In min_order  : minmun order
!   In max_order  : maximum order
!   In maxAutorized : maximum autorized order
!
! --------------------------------------------------------------------------------------------------
!
        ASSERT(min_order <= max_order)
!
        if (max_order > maxAutorized) then
            call utmess('F', 'HHO1_10', ni=2, vali=(/max_order, maxAutorized/))
        end if
!
        if (min_order < 0) then
            call utmess('F', 'HHO1_11', si=min_order)
        end if
!
    end subroutine
!
! -- member functions
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoBasisCellInit(this, hhoCell, type)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        class(HHO_basis_cell), intent(out)  :: this
        integer, optional, intent(in)       :: type
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   initialize hho basis for a cell
!   In hhoCell  : the current HHO cell
!   In type     : type of the basis
!   Out this    : HHO_basis_cell
!
! --------------------------------------------------------------------------------------------------
!
        integer :: idim, size_basis_scal, ipg, iret, jtab(1), ib, offset, size_face
        real(kind=8) :: axes(3, 3), length_box(3)
        real(kind=8), dimension(MSIZE_CELL_SCAL, MAX_QP_CELL) :: basisOrthoIpg
        type(HHO_basis_cell) :: hhoBasisIner
        type(HHO_Quadrature) :: hhoQuad
!
        call this%hhoMono%initialize(hhoCell%ndim, MAX_DEGREE_CELL)
!
        if (present(type)) then
            this%type = type
        end if
!
        if (this%type == BASIS_CARTESIAN) then
            axes = 0.d0
            do idim = 1, hhoCell%ndim
                axes(idim, idim) = 1.d0
            end do
        else
            axes = hhoCell%axes
        end if
        length_box = hhoLengthBoundingBoxCell(hhoCell, axes)
!
        this%rotmat = transpose(axes)
        this%ndim = hhoCell%ndim
        this%scaling_factor = 2.d0/length_box
        this%center = hhoCenterBoundingBoxCell(hhoCell, axes)
!
        do idim = 1, this%ndim
            this%rotmat(idim, :) = &
                this%rotmat(idim, :)*this%scaling_factor(idim)
        end do
!
        if (this%type == BASIS_ORTHO) then
!
            call hhoBasisIner%initialize(hhoCell, BASIS_INERTIAL)
            size_basis_scal = hhoBasisIner%BSSize(0, MAX_DEGREE_CELL)
!
            call tecach('NNO', 'PCHHOBS', 'L', iret, nval=1, itab=jtab)
!
            if (iret == 0) then
                this%coeff_shift(1) = 1
                do ib = 1, size_basis_scal
                    this%coeff_shift(ib+1) = this%coeff_shift(ib)+ib
                end do
                size_face = binomial(MAX_DEGREE_FACE+this%ndim-1, MAX_DEGREE_FACE)
                offset = hhoCell%nbfaces*(size_face*(size_face+1)/2)
                call readVector('PCHHOBS', maxval(this%coeff_shift)-1, &
                                this%coeff_mono, offset)
            else
!
! ------------ If you have this error - add the basis field as an input of you option
                call jevech('PCHHOBO', 'E', iret)

                basisOrthoIpg = 0.d0
!
                if (this%ndim > 1) then
!               Ortho-normalisation with modified Gramm-Schmidt
!
                    call hhoQuad%GetQuadCell(hhoCell, 2*MAX_DEGREE_CELL)
!
! --------------------- Initialize phi_i
!
                    do ipg = 1, hhoQuad%nbQuadPoints
                        call hhoBasisIner%BSEval(hhoQuad%points(1:3, ipg), &
                                                 0, MAX_DEGREE_CELL, basisOrthoIpg(1, ipg))
                    end do
                end if
!
                call orthonormalization(hhoQuad, basisOrthoIpg, size_basis_scal, this%ndim, &
                                        hhoCell%measure, this%coeff_shift, this%coeff_mono)
!
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoBasisFaceInit(this, hhoFace, type)
!
        implicit none
!
        type(HHO_Face), intent(in)               :: hhoFace
        class(HHO_basis_face), intent(out)       :: this
        integer, optional, intent(in)            :: type
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   initialize hho basis for a face
!   In hhoFace     : the current HHO face
!   In type        : type of the basis (default: AUTO)
!   Out this       : HHO_basis_face
!
! --------------------------------------------------------------------------------------------------
!
        integer :: idim, size_basis_scal, ipg, iret, jtab(1), ib, offset, nb_coeff
        real(kind=8) :: axes(3, 2), length_box(2)
        real(kind=8), dimension(MSIZE_CELL_SCAL, MAX_QP_FACE) :: basisOrthoIpg
        type(HHO_basis_face) :: hhoBasisIner
        type(HHO_Quadrature) :: hhoQuad
!
        call this%hhoMono%initialize(hhoFace%ndim, MAX_DEGREE_FACE)
!
        if (present(type)) then
            this%type = type
        end if
!
        if (this%type == BASIS_CARTESIAN) then
            axes = hhoLocalBasisFace(hhoFace)
        else
            axes = hhoFace%axes
        end if
        length_box = hhoLengthBoundingBoxFace(hhoFace, axes)
!
        this%ndim = hhoFace%ndim
        this%scaling_factor = 2.d0/length_box
        this%rotmat = transpose(axes)
        this%center = hhoCenterBoundingBoxFace(hhoFace, axes)
!
        do idim = 1, this%ndim
            this%rotmat(idim, :) = &
                this%rotmat(idim, :)*this%scaling_factor(idim)
        end do
!
        if (this%type == BASIS_ORTHO) then
!
            call hhoBasisIner%initialize(hhoFace, BASIS_INERTIAL)
!
            size_basis_scal = hhoBasisIner%BSSize(0, MAX_DEGREE_FACE)
!
            call tecach('NNO', 'PCHHOBS', 'L', iret, nval=1, itab=jtab)
!
            if (iret == 0) then
                this%coeff_shift(1) = 1
                do ib = 1, size_basis_scal
                    this%coeff_shift(ib+1) = this%coeff_shift(ib)+ib
                end do
!
                nb_coeff = this%coeff_shift(size_basis_scal+1)-1
                offset = (hhoFace%face_loc-1)*nb_coeff
                call readVector('PCHHOBS', nb_coeff, this%coeff_mono, offset)
!
            else
!
! ------------ If you have this error - add the basis field as an input of your option
                call jevech('PCHHOBO', 'E', iret)
!
                basisOrthoIpg = 0.d0
!
                if (this%ndim > 1) then
!
                    call hhoQuad%GetQuadFace(hhoFace, 2*MAX_DEGREE_FACE)
!
! --------------------- Initialize phi_i
!
                    do ipg = 1, hhoQuad%nbQuadPoints
                        call hhoBasisIner%BSEval(hhoQuad%points(1:3, ipg), &
                                                 0, MAX_DEGREE_FACE, basisOrthoIpg(1, ipg))
                    end do
                end if
!
                call orthonormalization(hhoQuad, basisOrthoIpg, size_basis_scal, this%ndim, &
                                        hhoFace%measure, this%coeff_shift, this%coeff_mono)
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoBSCellSize(this, min_order, max_order) result(size_basis)
!
        implicit none
!
        class(HHO_basis_cell), intent(in)      :: this
        integer, intent(in)                    :: min_order
        integer, intent(in)                    :: max_order
        integer                                :: size_basis
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the basis for the specified range [min_order, max_order]
!   In this                 : HHO Basis cell
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out size_basis          : size of the scalar basis
! --------------------------------------------------------------------------------------------------
!
        integer :: size_min, size_max
!
        call this%BSRange(min_order, max_order, size_min, size_max)
        size_basis = size_max-size_min+1
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoBSFaceSize(this, min_order, max_order) result(size_basis)
!
        implicit none
!
        class(HHO_basis_face), intent(in)       :: this
        integer, intent(in)                     :: min_order
        integer, intent(in)                     :: max_order
        integer                                 :: size_basis
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the basis for the specified range [min_order, max_order]
!   In this                 : HHO Basis Face
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out size_basis          : size of the scalar basis
! --------------------------------------------------------------------------------------------------
!
        integer :: size_min, size_max
!
        call this%BSRange(min_order, max_order, size_min, size_max)
        size_basis = size_max-size_min+1
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoBSCellRange(this, min_order, max_order, ifrom, ito)
!
        implicit none
!
        class(HHO_basis_cell), intent(in)       :: this
        integer, intent(in)                     :: min_order
        integer, intent(in)                     :: max_order
        integer, intent(out)                    :: ifrom
        integer, intent(out)                    :: ito
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the basis for the specified range [min_order, max_order]
!   In this                 : HHO Basis cell
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out ifrom               : first index of the monomials
!   Out ito                 : last index of the monomials
! --------------------------------------------------------------------------------------------------
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
        if (min_order == 0) then
            ifrom = 1
        else
            ifrom = binomial(min_order-1+this%hhoMono%ndim, min_order-1)+1
        end if
!
        ito = binomial(max_order+this%hhoMono%ndim, max_order)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoBSFaceRange(this, min_order, max_order, ifrom, ito)
!
        implicit none
!
        class(HHO_basis_face), intent(in)       :: this
        integer, intent(in)                     :: min_order
        integer, intent(in)                     :: max_order
        integer, intent(out)                    :: ifrom
        integer, intent(out)                    :: ito
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the basis for the specified range [min_order, max_order]
!   In this                 : HHO BS Face
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out ifrom               : first index of the monomials
!   Out ito                 : last index of the monomials
! --------------------------------------------------------------------------------------------------
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
        if (min_order == 0) then
            ifrom = 1
        else
            ifrom = binomial(min_order-1+this%hhoMono%ndim, min_order-1)+1
        end if
!
        ito = binomial(max_order+this%hhoMono%ndim, max_order)
!
    end subroutine
!
    function hhoBVCellSize(this, min_order, max_order) result(size_basis)
!
!===================================================================================================
!
!===================================================================================================
!
        implicit none
!
        class(HHO_basis_cell), intent(in)      :: this
        integer, intent(in)                    :: min_order
        integer, intent(in)                    :: max_order
        integer                                :: size_basis
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the basis for the specified range [min_order, max_order]
!   In this                 : HHO Basis cell
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out size_basis          : size of the scalar basis
! --------------------------------------------------------------------------------------------------
!
        integer :: size_min, size_max
!
        call this%BSRange(min_order, max_order, size_min, size_max)
        size_basis = this%hhoMono%ndim*(size_max-size_min+1)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoBVFaceSize(this, min_order, max_order) result(size_basis)
!
        implicit none
!
        class(HHO_basis_face), intent(in)       :: this
        integer, intent(in)                     :: min_order
        integer, intent(in)                     :: max_order
        integer                                 :: size_basis
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the basis for the specified range [min_order, max_order]
!   In this                 : HHO Basis Face
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out size_basis          : size of the scalar basis
! --------------------------------------------------------------------------------------------------
!
        integer :: size_min, size_max
!
        call this%BSRange(min_order, max_order, size_min, size_max)
        size_basis = (this%hhoMono%ndim+1)*(size_max-size_min+1)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoBVCellRange(this, min_order, max_order, ifrom, ito)
!
        implicit none
!
        class(HHO_basis_cell), intent(in)       :: this
        integer, intent(in)                     :: min_order
        integer, intent(in)                     :: max_order
        integer, intent(out)                    :: ifrom
        integer, intent(out)                    :: ito
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the basis for the specified range [min_order, max_order]
!   In this                 : HHO Basis cell
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out ifrom               : first index of the monomials
!   Out ito                 : last index of the monomials
! --------------------------------------------------------------------------------------------------
!
        integer :: ndim
!
        ndim = this%hhoMono%ndim
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
        if (min_order == 0) then
            ifrom = 1
        else
            ifrom = ndim*binomial(min_order-1+ndim, min_order-1)+1
        end if
!
        ito = ndim*binomial(max_order+ndim, max_order)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoBVFaceRange(this, min_order, max_order, ifrom, ito)
!
        implicit none
!
        class(HHO_basis_face), intent(in)       :: this
        integer, intent(in)                     :: min_order
        integer, intent(in)                     :: max_order
        integer, intent(out)                    :: ifrom
        integer, intent(out)                    :: ito
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the basis for the specified range [min_order, max_order]
!   In this                 : HHO BS Face
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out ifrom               : first index of the monomials
!   Out ito                 : last index of the monomials
! --------------------------------------------------------------------------------------------------
!
        integer :: ndim
!
        ndim = this%hhoMono%ndim
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
        if (min_order == 0) then
            ifrom = 1
        else
            ifrom = (ndim+1)*binomial(min_order-1+ndim, min_order-1)+1
        end if
        ito = (ndim+1)*binomial(max_order+ndim, max_order)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoBMCellSize(this, min_order, max_order) result(size_basis)
!
        implicit none
!
        class(HHO_basis_cell), intent(in)      :: this
        integer, intent(in)                    :: min_order
        integer, intent(in)                    :: max_order
        integer                                :: size_basis
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the matrix basis for the specified range [min_order, max_order]
!   In this                 : HHO Basis cell
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out size_basis          : size of the matrix basis
! --------------------------------------------------------------------------------------------------
!
        integer :: size_min, size_max, ndim2
!
        ndim2 = this%hhoMono%ndim*this%hhoMono%ndim
!
        call this%BSRange(min_order, max_order, size_min, size_max)
        size_basis = ndim2*(size_max-size_min+1)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoBMCellRange(this, min_order, max_order, ifrom, ito)
!
        implicit none
!
        class(HHO_basis_cell), intent(in)       :: this
        integer, intent(in)                     :: min_order
        integer, intent(in)                     :: max_order
        integer, intent(out)                    :: ifrom
        integer, intent(out)                    :: ito
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   size of the matrix basis for the specified range [min_order, max_order]
!   In this                 : HHO Basis cell
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out ifrom               : first index of the monomials
!   Out ito                 : last index of the monomials
! --------------------------------------------------------------------------------------------------
!
        integer :: ndim, ndim2
!
        ndim = this%hhoMono%ndim
        ndim2 = ndim*ndim
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
        if (min_order == 0) then
            ifrom = 1
        else
            ifrom = ndim2*binomial(min_order-1+ndim, min_order-1)+1
        end if
        ito = ndim2*binomial(max_order+ndim, max_order)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoBSCellEval(this, point, min_order, max_order, basisScalEval)
!
        implicit none
!
        class(HHO_basis_cell), intent(inout)                    :: this
        real(kind=8), dimension(3), intent(in)                  :: point
        integer, intent(in)                                     :: min_order
        integer, intent(in)                                     :: max_order
        real(kind=8), dimension(MSIZE_CELL_SCAL), intent(out)   :: basisScalEval
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   evaluate hho basis scalar for a 3D cell
!   In this                 : HHO_basis_scalar_cell
!   In point                : point where evaluate
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out basisScalEval       : evaluation of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: peval
        integer :: imono, ifrom, ito, icoeff
        integer, dimension(3) :: power
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
        call this%BSRange(min_order, max_order, ifrom, ito)
! ----  scaled point
        peval = this%map_pt(point)
!
! ----- Eval monomials
        call this%hhoMono%eval(peval)
!
        basisScalEval = 0.d0
!
! ----- Loop on the monomial
!
        if (this%ndim == 1) then
            if (this%type == BASIS_ORTHO) then
                do imono = ifrom, ito
                    do icoeff = 1, this%coeff_shift(imono+1)-this%coeff_shift(imono)
                        power(1) = this%hhoMono%monomials(1, icoeff)
                        basisScalEval(imono) = basisScalEval(imono)+ &
                                               this%coeff_mono(this%coeff_shift(imono)-1+icoeff)* &
                                               this%hhoMono%monoEval(1, power(1))
                    end do
                end do
            else
                do imono = ifrom, ito
                    power(1) = this%hhoMono%monomials(1, imono)
                    basisScalEval(imono) = this%hhoMono%monoEval(1, power(1))
                end do
            end if
        else if (this%ndim == 2) then
            if (this%type == BASIS_ORTHO) then
                do imono = ifrom, ito
                    do icoeff = 1, this%coeff_shift(imono+1)-this%coeff_shift(imono)
                        power(1:2) = this%hhoMono%monomials(1:2, icoeff)
                        basisScalEval(imono) = basisScalEval(imono)+ &
                                               this%coeff_mono(this%coeff_shift(imono)-1+icoeff)* &
                                               this%hhoMono%monoEval(1, power(1))* &
                                               this%hhoMono%monoEval(2, power(2))
                    end do
                end do
            else
                do imono = ifrom, ito
                    power(1:2) = this%hhoMono%monomials(1:2, imono)
                    basisScalEval(imono) = this%hhoMono%monoEval(1, power(1))* &
                                           this%hhoMono%monoEval(2, power(2))
                end do
            end if
        else if (this%ndim == 3) then
            if (this%type == BASIS_ORTHO) then
                do imono = ifrom, ito
                    do icoeff = 1, this%coeff_shift(imono+1)-this%coeff_shift(imono)
                        power(1:3) = this%hhoMono%monomials(1:3, icoeff)
                        basisScalEval(imono) = basisScalEval(imono)+ &
                                               this%coeff_mono(this%coeff_shift(imono)-1+icoeff)* &
                                               this%hhoMono%monoEval(1, power(1))* &
                                               this%hhoMono%monoEval(2, power(2))* &
                                               this%hhoMono%monoEval(3, power(3))
                    end do
                end do
            else
                do imono = ifrom, ito
                    power(1:3) = this%hhoMono%monomials(1:3, imono)
                    basisScalEval(imono) = this%hhoMono%monoEval(1, power(1))* &
                                           this%hhoMono%monoEval(2, power(2))* &
                                           this%hhoMono%monoEval(3, power(3))

                end do
            end if
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
    subroutine hhoBSCellGradEv(this, point, min_order, max_order, BSGradEval)
!
        implicit none
!
        class(HHO_basis_cell), intent(inout)                    :: this
        real(kind=8), dimension(3), intent(in)                  :: point
        integer, intent(in)                                     :: min_order
        integer, intent(in)                                     :: max_order
        real(kind=8), dimension(3, MSIZE_CELL_SCAL), intent(out):: BSGradEval
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   evaluate hho basis scalar for a 3D cell
!   In this                 : HHO_basis_scalar_cell
!   In point                : point where evaluate
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out BSGradEval   : evaluation of the gradient of the scalar basis
!
! --------------------------------------------------------------------------------------------------
!
        integer :: ifrom, ito, icoeff
        real(kind=8), dimension(3) :: peval, func, dfunc, grad
        integer :: ind, imono
        integer, dimension(3) :: power
        real(kind=8) :: invrotmat(3, 3)
        real(kind=8), dimension(3, MSIZE_CELL_SCAL) :: Grad_mono
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
        call this%BSRange(min_order, max_order, ifrom, ito)
!
        invrotmat = transpose(this%rotmat)
! ----  scaled point
        peval = this%map_pt(point)
!
! ----- Eval monomials
        call this%hhoMono%eval(peval)
!
        BSGradEval = 0.d0
        Grad_mono = 0.d0
        grad = 0.d0
!
! ----- Loop on the monomial
        ind = 0
        if (this%ndim == 1) then
            if (this%type == BASIS_ORTHO) then
                do imono = 1, ifrom-1
                    power(1) = this%hhoMono%monomials(1, imono)
!
                    func(1) = this%hhoMono%monoEval(1, power(1))
!
                    if (power(1) == 0) then
                        dfunc(1) = 0.d0
                    else
                        dfunc(1) = power(1)*this%hhoMono%monoEval(1, power(1)-1)* &
                                   this%scaling_factor(1)
                    end if
!
                    Grad_mono(1, imono) = dfunc(1)
                end do
            end if
!
            do imono = ifrom, ito
                ind = ind+1
                power(1) = this%hhoMono%monomials(1, imono)
!
                func(1) = this%hhoMono%monoEval(1, power(1))
!
                if (power(1) == 0) then
                    dfunc(1) = 0.d0
                else
                    dfunc(1) = power(1)*this%hhoMono%monoEval(1, power(1)-1)* &
                               this%scaling_factor(1)
                end if
!
                BSGradEval(1, ind) = dfunc(1)
!
                if (this%type == BASIS_ORTHO) then
                    Grad_mono(1, ind) = BSGradEval(1, ind)
                    BSGradEval(1, ind) = 0.d0
                    do icoeff = 1, this%coeff_shift(imono+1)-this%coeff_shift(imono)
                        BSGradEval(1, ind) = BSGradEval(1, ind)+ &
                                             this%coeff_mono(this%coeff_shift(imono)-1+icoeff)* &
                                             Grad_mono(1, icoeff)
                    end do
                end if
            end do
!
        else if (this%ndim == 2) then
            if (this%type == BASIS_ORTHO) then
                do imono = 1, ifrom-1
                    power(1:2) = this%hhoMono%monomials(1:2, imono)
!
                    func(1:2) = (/this%hhoMono%monoEval(1, power(1)), &
                                  this%hhoMono%monoEval(2, power(2))/)
!
                    if (power(1) == 0) then
                        dfunc(1) = 0.d0
                    else
                        dfunc(1) = power(1)*this%hhoMono%monoEval(1, power(1)-1)
                    end if
                    if (power(2) == 0) then
                        dfunc(2) = 0.d0
                    else
                        dfunc(2) = power(2)*this%hhoMono%monoEval(2, power(2)-1)
                    end if
!
                    grad(1) = dfunc(1)*func(2)
                    grad(2) = func(1)*dfunc(2)
                    Grad_mono(1:2, imono) = matmul(invrotmat(1:2, 1:2), grad(1:2))
                end do
            end if
!
            do imono = ifrom, ito
                ind = ind+1
                power(1:2) = this%hhoMono%monomials(1:2, imono)
!
                func(1:2) = (/this%hhoMono%monoEval(1, power(1)), &
                              this%hhoMono%monoEval(2, power(2))/)
!
                if (power(1) == 0) then
                    dfunc(1) = 0.d0
                else
                    dfunc(1) = power(1)*this%hhoMono%monoEval(1, power(1)-1)
                end if
                if (power(2) == 0) then
                    dfunc(2) = 0.d0
                else
                    dfunc(2) = power(2)*this%hhoMono%monoEval(2, power(2)-1)
                end if
!
                grad(1) = dfunc(1)*func(2)
                grad(2) = func(1)*dfunc(2)
                BSGradEval(1:2, ind) = matmul(invrotmat(1:2, 1:2), grad(1:2))
!
                if (this%type == BASIS_ORTHO) then
                    Grad_mono(1:2, imono) = BSGradEval(1:2, ind)
                    BSGradEval(1:2, ind) = 0.d0
                    do icoeff = 1, this%coeff_shift(imono+1)-this%coeff_shift(imono)
                        BSGradEval(1:2, ind) = BSGradEval(1:2, ind)+ &
                                               this%coeff_mono(this%coeff_shift(imono)-1+icoeff)* &
                                               Grad_mono(1:2, icoeff)
                    end do
                end if
            end do
!
        else if (this%ndim == 3) then
            if (this%type == BASIS_ORTHO) then
                do imono = 1, ifrom-1

                    power(1:3) = this%hhoMono%monomials(1:3, imono)
!
                    func(1:3) = (/this%hhoMono%monoEval(1, power(1)), &
                                  this%hhoMono%monoEval(2, power(2)), &
                                  this%hhoMono%monoEval(3, power(3))/)
!
                    if (power(1) == 0) then
                        dfunc(1) = 0.d0
                    else
                        dfunc(1) = power(1)*this%hhoMono%monoEval(1, power(1)-1)
                    end if
                    if (power(2) == 0) then
                        dfunc(2) = 0.d0
                    else
                        dfunc(2) = power(2)*this%hhoMono%monoEval(2, power(2)-1)
                    end if
                    if (power(3) == 0) then
                        dfunc(3) = 0.d0
                    else
                        dfunc(3) = power(3)*this%hhoMono%monoEval(3, power(3)-1)
                    end if
!
                    grad(1) = dfunc(1)*func(2)*func(3)
                    grad(2) = func(1)*dfunc(2)*func(3)
                    grad(3) = func(1)*func(2)*dfunc(3)
                    Grad_mono(1:3, imono) = matmul(invrotmat, grad)
                end do
            end if
!
            do imono = ifrom, ito
                ind = ind+1
                power(1:3) = this%hhoMono%monomials(1:3, imono)
!
                func(1:3) = (/this%hhoMono%monoEval(1, power(1)), &
                              this%hhoMono%monoEval(2, power(2)), &
                              this%hhoMono%monoEval(3, power(3))/)
!
                if (power(1) == 0) then
                    dfunc(1) = 0.d0
                else
                    dfunc(1) = power(1)*this%hhoMono%monoEval(1, power(1)-1)
                end if
                if (power(2) == 0) then
                    dfunc(2) = 0.d0
                else
                    dfunc(2) = power(2)*this%hhoMono%monoEval(2, power(2)-1)
                end if
                if (power(3) == 0) then
                    dfunc(3) = 0.d0
                else
                    dfunc(3) = power(3)*this%hhoMono%monoEval(3, power(3)-1)
                end if
!
                grad(1) = dfunc(1)*func(2)*func(3)
                grad(2) = func(1)*dfunc(2)*func(3)
                grad(3) = func(1)*func(2)*dfunc(3)
                BSGradEval(1:3, ind) = matmul(invrotmat, grad)

                if (this%type == BASIS_ORTHO) then
                    Grad_mono(1:3, imono) = BSGradEval(1:3, ind)
                    BSGradEval(1:3, ind) = 0.d0
                    do icoeff = 1, this%coeff_shift(imono+1)-this%coeff_shift(imono)
                        BSGradEval(1:3, ind) = BSGradEval(1:3, ind)+ &
                                               this%coeff_mono(this%coeff_shift(imono)-1+icoeff)* &
                                               Grad_mono(1:3, icoeff)
                    end do
                end if
            end do
!
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
    subroutine hhoBVCellSymGdEv(this, point, min_order, max_order, BVSymGradEval)
!
        implicit none
!
        class(HHO_basis_cell), intent(inout)                    :: this
        real(kind=8), dimension(3), intent(in)                  :: point
        integer, intent(in)                                     :: min_order
        integer, intent(in)                                     :: max_order
        real(kind=8), dimension(6, MSIZE_CELL_VEC), intent(out) :: BVSymGradEval
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   evaluate hho gradient basis vectoriel for a 3D cell
!   In this                 : HHO_basis_cell
!   In point                : point where evaluate
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out BVSymGradEval       : evaluation of the symmetric gradient of the vectoriel basis
!
!   Be carefull the format of symetric matrix M is a vector of six components
!   (M11, M22, M33, Rac2*M12, Rac2*M13, Rac2*M23)
!
! --------------------------------------------------------------------------------------------------
!
        integer :: size_basis_scal, ind, imono
        real(kind=8), parameter :: rac2sur2 = sqrt(2.d0)/2.d0
        real(kind=8), dimension(3, MSIZE_CELL_SCAL) :: BSGradEval
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
        size_basis_scal = this%BSSize(min_order, max_order)
!
! ----- Eval scalar gradient
!
        call this%BSEvalGrad(point, min_order, max_order, BSGradEval)
!
! ----- Loop on the monomial
        BVSymGradEval = 0.d0
        ind = 0
        if (this%ndim == 2) then
! -------- dir = 1
            do imono = 1, size_basis_scal
                ind = ind+1
!
                BVSymGradEval(1, ind) = BSGradEval(1, imono)
                BVSymGradEval(4, ind) = BSGradEval(2, imono)*rac2sur2
            end do
! -------- dir = 2
            do imono = 1, size_basis_scal
                ind = ind+1
!
                BVSymGradEval(2, ind) = BSGradEval(2, imono)
                BVSymGradEval(4, ind) = BSGradEval(1, imono)*rac2sur2
            end do
        else if (this%ndim == 3) then
! -------- dir = 1
            do imono = 1, size_basis_scal
                ind = ind+1
!
                BVSymGradEval(1, ind) = BSGradEval(1, imono)
                BVSymGradEval(4, ind) = BSGradEval(2, imono)*rac2sur2
                BVSymGradEval(5, ind) = BSGradEval(3, imono)*rac2sur2
            end do
! -------- dir = 2
            do imono = 1, size_basis_scal
                ind = ind+1
!
                BVSymGradEval(2, ind) = BSGradEval(2, imono)
                BVSymGradEval(4, ind) = BSGradEval(1, imono)*rac2sur2
                BVSymGradEval(6, ind) = BSGradEval(3, imono)*rac2sur2
            end do
! -------- dir = 3
            do imono = 1, size_basis_scal
                ind = ind+1
!
                BVSymGradEval(3, ind) = BSGradEval(3, imono)
                BVSymGradEval(5, ind) = BSGradEval(1, imono)*rac2sur2
                BVSymGradEval(6, ind) = BSGradEval(2, imono)*rac2sur2
            end do
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
    subroutine hhoBSFaceEval(this, point, min_order, max_order, basisScalEval)
!
        implicit none
!
        class(HHO_basis_face), intent(inout)                    :: this
        real(kind=8), dimension(3), intent(in)                  :: point
        integer, intent(in)                                     :: min_order
        integer, intent(in)                                     :: max_order
        real(kind=8), dimension(MSIZE_FACE_SCAL), intent(out)   :: basisScalEval
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   evaluate hho basis scalar function for a face
!   In this                 : HHO_basis_face
!   In point                : point where evaluate
!   In min_order            : minimum order
!   In max_order            : maximum order
!   Out basisScalEval       : evaluation of the scalar basis function
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(2) :: peval
        integer :: imono, ifrom, ito, icoeff
        integer, dimension(2) :: power
!
! ---- Check the order
        call check_order(min_order, max_order, this%hhoMono%maxOrder())
!
! ----  scaled point
        peval = this%map_pt(point)
!
! ----- Eval monomials
        call this%hhoMono%eval([peval(1), peval(2), 0.d0])
!
        basisScalEval = 0.d0
        call this%BSRange(min_order, max_order, ifrom, ito)
!
! ----- Loop on the monomial
        if (this%ndim == 0) then
            basisScalEval(1) = 1.d0
        else if (this%ndim == 1) then
            if (this%type == BASIS_ORTHO) then
                do imono = ifrom, ito
                    do icoeff = 1, this%coeff_shift(imono+1)-this%coeff_shift(imono)
                        power(1) = this%hhoMono%monomials(1, icoeff)
                        basisScalEval(imono) = basisScalEval(imono)+ &
                                               this%coeff_mono(this%coeff_shift(imono)-1+icoeff)* &
                                               this%hhoMono%monoEval(1, power(1))
                    end do
                end do
            else
                do imono = ifrom, ito
                    power(1) = this%hhoMono%monomials(1, imono)
                    basisScalEval(imono) = this%hhoMono%monoEval(1, power(1))
                end do
            end if
        else if (this%ndim == 2) then
            if (this%type == BASIS_ORTHO) then
                do imono = ifrom, ito
                    do icoeff = 1, this%coeff_shift(imono+1)-this%coeff_shift(imono)
                        power(1:2) = this%hhoMono%monomials(1:2, icoeff)
                        basisScalEval(imono) = basisScalEval(imono)+ &
                                               this%coeff_mono(this%coeff_shift(imono)-1+icoeff)* &
                                               this%hhoMono%monoEval(1, power(1))* &
                                               this%hhoMono%monoEval(2, power(2))
                    end do
                end do
            else
                do imono = ifrom, ito
                    power(1:2) = this%hhoMono%monomials(1:2, imono)
                    basisScalEval(imono) = this%hhoMono%monoEval(1, power(1))* &
                                           this%hhoMono%monoEval(2, power(2))
                end do
            end if
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
    function map_pt_face(this, point) result(proj)
!
        implicit none
!
        class(HHO_basis_face), intent(in)          :: this
        real(kind=8), dimension(3), intent(in)     :: point
        real(kind=8), dimension(2)                 :: proj
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   map a point in global coordinates to local coordinates
!   In point                : point to project
!   Out proj                : projection
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: ep
!
        ep(1:3) = point(1:3)-this%center(1:3)
!
        proj = matmul(this%rotmat, ep)
        ! if(proj(1) < -1.1d0 .or. proj(1) > 1.1d0) print*, "errorf0: ", proj
        ! if(proj(2) < -1.1d0 .or. proj(2) > 1.1d0) print*, "errorf1: ", proj
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function map_pt_cell(this, point) result(proj)
!
        implicit none
!
        class(HHO_basis_cell), intent(in)          :: this
        real(kind=8), dimension(3), intent(in)     :: point
        real(kind=8), dimension(3)                 :: proj
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   map a point in global coordinates to local coordinates
!   In point                : point to project
!   Out proj                : projection
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: ep
!
        ep(1:3) = point(1:3)-this%center(1:3)
!
        proj = matmul(this%rotmat, ep)
        ! if(proj(1) < -1.1d0 .or. proj(1) > 1.1d0) print*, "errorc0: ", proj
        ! if(proj(2) < -1.1d0 .or. proj(2) > 1.1d0) print*, "errorc1: ", proj
        ! if(proj(3) < -1.1d0 .or. proj(3) > 1.1d0) print*, "errorc2: ", proj
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine orthonormalization(hhoQuad, basisIpg, nb_basis, ndim, measure, &
                                  coeff_shift, coeff_mono)
!
        implicit none
!
        type(HHO_Quadrature), intent(in)                    :: hhoQuad
        real(kind=8), intent(inout)                         :: basisIpg(MSIZE_CELL_SCAL, *)
        integer, intent(in)                                 :: nb_basis, ndim
        real(kind=8), intent(in)                            :: measure
        integer, intent(out)                                :: coeff_shift(*)
        real(kind=8), intent(out)                           :: coeff_mono(*)
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   Modified Gramm-Schmidt algorithm to orthonormalize basis functions
!   In hhoQuad              : Quadrature to compute scalar product
!   IO basisIpg             : basis to orthonormalize (at quadrature point)
!   In nb_basis             : number of basis function
!   In ndim                 : dimension of the basis
!   Out coeff_mono          : coefficient of the orthonormal basis function
!   Out coeff_shift         : shifting for each basis function
!
! --------------------------------------------------------------------------------------------------
!
        integer, parameter :: nb_ortho = 2
        integer :: i, j, k, ipg, npg, i_ortho
        real(kind=8) :: rp(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL), rc(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
        real(kind=8) :: ra(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
        real(kind=8) :: ri(MSIZE_CELL_SCAL)
!
        ASSERT(nb_basis <= MSIZE_CELL_SCAL)
!
        if (ndim == 1) then
            ! Coeffficient of Legendre basis
            coeff_shift(1:4) = [1, 2, 4, 7]
            coeff_mono(1:6) = &
                ! ordre 0
                (/1.d0, &
                  ! ordre 1
                  0.d0, 1.d0, &
                  ! ordre 2
                  -0.5d0, 0.d0, 1.5d0/)
!
            coeff_mono(1) = coeff_mono(1)/sqrt(measure)
            coeff_mono(2:3) = coeff_mono(2:3)*sqrt(3.d0/measure)
            coeff_mono(4:6) = coeff_mono(4:6)*sqrt(5.d0/measure)
!
        else
!           Ortho-normalisation with modified Gramm-Schmidt
!
            coeff_shift(1) = 1
            do i = 2, nb_basis+1
                coeff_shift(i) = coeff_shift(i-1)+i-1
            end do
!
            npg = hhoQuad%nbQuadPoints
!
            rp = 0.d0
            do i_ortho = 1, nb_ortho
                rc = 0.d0
                ra = 0.d0

                do i = 1, nb_basis
!
                    ri = 0.d0
                    do j = 1, i-1
!
! --------------------- Compute r_ij = (phi_i, phi_j)_T
!
                        do ipg = 1, npg
                            ri(j) = ri(j)+hhoQuad%weights(ipg)* &
                                    basisIpg(i, ipg)*basisIpg(j, ipg)
                        end do
!
! --------------------- Update phi_i
!
                        do ipg = 1, npg
                            basisIpg(i, ipg) = basisIpg(i, ipg)-ri(j)*basisIpg(j, ipg)
                        end do
                    end do
!
! ------------------ Compute normalization
!
                    do ipg = 1, npg
                        ri(i) = ri(i)+hhoQuad%weights(ipg)*basisIpg(i, ipg)*basisIpg(i, ipg)
                    end do
!
! ------------------ Rescale coefficient
!
                    ri(i) = 1.d0/sqrt(ri(i))
                    ri(1:i-1) = ri(1:i-1)*ri(i)
                    basisIpg(i, 1:npg) = basisIpg(i, 1:npg)*ri(i)
!
! ------------- Compute coefficient of orthogonal basis function
!
                    rc(i, i) = ri(i)
                    do j = 1, i-1
                        do k = 1, j
                            rc(i, k) = rc(i, k)-ri(j)*rc(j, k)
                        end do
                    end do

                    if (i_ortho == 1) then
                        ra(i, 1:i) = rc(i, 1:i)
                    else
                        do j = 1, i
                            do k = j, i
                                ra(i, j) = ra(i, j)+rc(i, k)*rp(k, j)
                            end do
                        end do
                    end if
                end do
                rp(1:nb_basis, 1:nb_basis) = ra(1:nb_basis, 1:nb_basis)
            end do

            do i = 1, nb_basis
                coeff_mono(coeff_shift(i)+i-1) = ra(i, i)
                do j = 1, i-1
                    coeff_mono(coeff_shift(i)+j-1) = ra(i, j)
                end do
            end do
!
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoBasisCellType(this) result(basis_type)
!
        implicit none
!
        class(HHO_basis_cell), intent(in)       :: this
        aster_logical                           :: basis_type
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   indicating the type of cell basis function used
!   Out basis_type  : booleen that indicates true if orthonormal
!
! --------------------------------------------------------------------------------------------------
!
        basis_type = this%type == BASIS_ORTHO
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoBasisFaceType(this) result(basis_type)
!
        implicit none
!
        class(HHO_basis_face), intent(in)       :: this
        aster_logical                           :: basis_type
!
! --------------------------------------------------------------------------------------------------
!   HHO - basis functions
!
!   indicating the type of face basis function used
!   Out basis_type  : booleen that indicates true if orthonormal
!
! --------------------------------------------------------------------------------------------------
!
        basis_type = this%type == BASIS_ORTHO
!
    end function
!
end module
