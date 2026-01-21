! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
module coupling_computation_module
!
    use coupling_type
    use FE_Basis_module
    use FE_mass_module
    use FE_quadrature_module
    use FE_topo_module
    use HHO_basis_module
    use HHO_massmat_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/coupling_penalisation_module.h"
#include "asterfort/HHO_size_module.h"
#include "blas/dsyr.h"
#include "FE_basis_module.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling methods - computation
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: cplCmpPenaFEHHOMatScal, cplCmpPenaFEHHOMatVec
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplCmpPenaFEHHOMatScal(feFace, hhoFace, hhoData, FEQuad, hhoQuad, coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: feFace
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Data), intent(in) :: hhoData
        type(FE_Quadrature), intent(in) :: FEQuad
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8) :: coef_pena
        type(HHO_matrix), intent(out) :: lhs
!
!===================================================================================================
!    FEM/HHO - Compute matrix: gamma/hF*(uh-uF, vh-vF)
!              unknowns (uh, uF)
!===================================================================================================
!
        type(HHO_basis_face) :: hhoBasisFace
        type(FE_Basis) :: FEBasis
        real(kind=8), dimension(MAX_BS) :: BSEval
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: basisScalEval
        real(kind=8), dimension(MAX_BS+MSIZE_FACE_SCAL) :: BS
        real(kind=8) :: coeff
        integer(kind=8) :: fbs, nbDoFsFE, nbDoFs, ipg
        blas_int, parameter :: b_one = to_blas_int(1)
        blas_int :: b_lda, b_n

!
! -- init face basis
        call hhoBasisFace%initialize(hhoFace)
        call FEBasis%initFace(feFace)
!
! -- number of dofs
        call hhoTherFaceDofs(hhoFace, hhoData, fbs)
        nbDoFsFE = FEBasis%size
        nbDoFs = nbDoFsFE+fbs
!
        call lhs%initialize(nbDoFs, nbDoFs, 0.d0)
        b_n = to_blas_int(lhs%nrows)
        b_lda = to_blas_int(lhs%max_nrows)
!
! -- Loop on quadrature point
        do ipg = 1, FEQuad%nbQuadPoints
! ----- Eval face basis function at the quadrature point
! -------- Compute uh
            BSEval = FEBasis%func(FEQuad%points_param(1:3, ipg))
! -------- Compute uF
            call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%face_degree(), &
                                     basisScalEval)
! -------- Assemble [uh, -uF]
            BS(1:nbDoFsFE) = BSEval(1:nbDoFsFE)
            BS(nbDoFsFE+1:nbDoFs) = -basisScalEval(1:fbs)
!
! --------  Eval

            coeff = FEQuad%weights(ipg)*coef_pena/hhoFace%diameter
            call dsyr('U', b_n, coeff, BS, b_one, &
                      lhs%m, b_lda)
        end do
!
        call lhs%copySymU()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpPenaFEHHOMatVec(feFace, hhoFace, hhoData, FEQuad, hhoQuad, coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: feFace
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Data), intent(in) :: hhoData
        type(FE_Quadrature), intent(in) :: FEQuad
        type(HHO_Quadrature), intent(in) :: hhoQuad
        real(kind=8) :: coef_pena
        type(HHO_matrix), intent(out) :: lhs
!
!===================================================================================================
!    FEM/HHO - Compute matrix: gamma/hF*(uh-uF, vh-vF)
!              unknowns (uh, uF)
!===================================================================================================
!
        type(HHO_matrix) :: lhs_scal
        integer(kind=8) :: ndim, ndofs, fbs, fbs_scal, nFE_scal, nFE
        integer(kind=8) :: idim, i, j, deca_row, deca_col

        call cplCmpPenaFEHHOMatScal(feFace, hhoFace, hhoData, FEQuad, hhoQuad, &
                                    coef_pena, lhs_scal)
!
        ndim = hhoFace%ndim+1
        ndofs = ndim*lhs_scal%nrows
        call hhoMecaFaceDofs(hhoFace, hhoData, fbs)
        fbs_scal = fbs/ndim
        nFE_scal = lhs_scal%nrows-fbs_scal
        nFE = nFE_scal*ndim
!
        call lhs%initialize(ndofs, ndofs, 0.d0)
!
        do idim = 1, ndim
            deca_col = idim
            do j = 1, nFE_scal
                deca_row = idim
                do i = 1, j
                    lhs%m(deca_row, deca_col) = lhs_scal%m(i, j)
                    deca_row = deca_row+ndim
                end do
                deca_col = deca_col+ndim
            end do
            deca_col = nFE+(idim-1)*fbs_scal+1
            do j = nFE_scal+1, lhs_scal%nrows
                deca_row = idim
                do i = 1, nFE_scal
                    lhs%m(deca_row, deca_col) = lhs_scal%m(i, j)
                    deca_row = deca_row+ndim
                end do
                deca_col = deca_col+1
            end do
            deca_col = nFE+(idim-1)*fbs_scal+1
            do j = nFE_scal+1, lhs_scal%nrows
                deca_row = nFE+(idim-1)*fbs_scal+1
                do i = nFE_scal+1, j
                    lhs%m(deca_row, deca_col) = lhs_scal%m(i, j)
                    deca_row = deca_row+1
                end do
                deca_col = deca_col+1
            end do
        end do
!
        call lhs%copySymU()
!
        call lhs_scal%free()
!
    end subroutine
!===================================================================================================
!
!===================================================================================================
!
end module
