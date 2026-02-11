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
    use HHO_utils_module, only: hhoCopySymPartMat
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/coupling_type.h"
#include "asterfort/HHO_size_module.h"
#include "blas/dger.h"
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
    public :: cplCmpPenaFEMHHOMatScal, cplCmpPenaFEMHHOMatVec
    public :: cplCmpPenaFEMFEMMatScal, cplCmpPenaFEMFEMMatVec
    public :: cplCmpLagrFEMFEMMatScal, cplCmpLagrFEMFEMMatVec

!
contains
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplCmpPenaFEMHHOMatScal(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                       coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        real(kind=8), intent(in) :: coef_pena
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
        call hhoBasisFace%initialize(hhoFaceMa)
        call FEBasis%initFace(FEFaceSl)
!
! -- number of dofs
        call hhoTherFaceDofs(hhoFaceMa, hhoDataMa, fbs)
        nbDoFsFE = FEBasis%size
        nbDoFs = nbDoFsFE+fbs
!
        call lhs%initialize(nbDoFs, nbDoFs, 0.d0)
        b_n = to_blas_int(lhs%nrows)
        b_lda = to_blas_int(lhs%max_nrows)
!
! -- Loop on quadrature point
        do ipg = 1, FEQuadSl%nbQuadPoints
! ----- Eval face basis function at the quadrature point
! -------- Compute uh
            BSEval = FEBasis%func(FEQuadSl%points_param(1:3, ipg))
! -------- Compute uF
            call hhoBasisFace%BSEval(hhoQuadMa%points(1:3, ipg), 0, hhoDataMa%face_degree(), &
                                     basisScalEval)
! -------- Assemble [uh, -uF]
            BS(1:nbDoFsFE) = BSEval(1:nbDoFsFE)
            BS(nbDoFsFE+1:nbDoFs) = -basisScalEval(1:fbs)
!
! --------  Eval

            coeff = FEQuadSl%weights(ipg)*coef_pena/hhoFaceMa%diameter
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
    subroutine cplCmpPenaFEMHHOMatVec(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                      coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl
        type(HHO_Face), intent(in) :: hhoFaceMa
        type(HHO_Data), intent(in) :: hhoDataMa
        type(FE_Quadrature), intent(in) :: FEQuadSl
        type(HHO_Quadrature), intent(in) :: hhoQuadMa
        real(kind=8), intent(in) :: coef_pena
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

        call cplCmpPenaFEMHHOMatScal(FEFaceSl, hhoFaceMa, hhoDataMa, FEQuadSl, hhoQuadMa, &
                                     coef_pena, lhs_scal)
!
        ndim = hhoFaceMa%ndim+1
        ndofs = ndim*lhs_scal%nrows
        call hhoMecaFaceDofs(hhoFaceMa, hhoDataMa, fbs)
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
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplCmpPenaFEMFEMMatScal(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                                       coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        real(kind=8), intent(in) :: coef_pena
        real(kind=8), intent(out) :: lhs(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    FEM/FEM - Compute matrix: gamma/hF*(uh1-uh2, vh1-vh2)
!              unknowns (uh1, uh2)
!===================================================================================================
!
        type(FE_Basis) :: FEBasisSl, FEBasisMa
        real(kind=8), dimension(MAX_BS) :: BSEvalSl, BSEvalMa
        real(kind=8), dimension(MAX_BS+MAX_BS) :: BS
        real(kind=8) :: coeff, hF
        integer(kind=8) :: nbDoFsFEMa, nbDoFsFESl, nbDoFs, ipg
        blas_int, parameter :: b_one = to_blas_int(1)
        blas_int :: b_lda, b_n
!
! -- init face basis
        call FEBasisSl%initFace(FEFaceSl)
        call FEBasisMa%initFace(FEFaceMa)
!
! -- number of dofs
        nbDoFsFESl = FEBasisSl%size
        nbDoFsFEMa = FEBasisMa%size
        nbDoFs = nbDoFsFESl+nbDoFsFEMa
!
        b_n = to_blas_int(nbDoFs)
        b_lda = to_blas_int(MSIZE_CPL_PENA_FEM)
!
        hF = FEFaceSl%diameter()
!
        lhs = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuadSl%nbQuadPoints
! ----- Eval face basis function at the quadrature point
! -------- Compute uh1
            BSEvalSl = FEBasisSl%func(FEQuadSl%points_param(1:3, ipg))
            BSEvalMa = FEBasisMa%func(FEQuadMa%points_param(1:3, ipg))
!
! -------- Assemble [uh1, -uh2]
            BS(1:nbDoFsFESl) = BSEvalSl(1:nbDoFsFESl)
            BS(nbDoFsFESl+1:nbDoFs) = -BSEvalMa(1:nbDoFsFEMa)
!
! --------  Eval

            coeff = FEQuadSl%weights(ipg)*coef_pena/hF
            call dsyr('U', b_n, coeff, BS, b_one, lhs, b_lda)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpPenaFEMFEMMatVec(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, coef_pena, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        real(kind=8), intent(in) :: coef_pena
        real(kind=8), intent(out) :: lhs(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    FEM/FEM - Compute matrix: gamma/hF*(uh1-uh2, vh1-vh2)
!              unknowns (uh1, uh2)
!===================================================================================================
!
        real(kind=8) :: lhs_scal(MSIZE_CPL_PENA_FEM, MSIZE_CPL_PENA_FEM)
        integer(kind=8) :: ndim, ndofsSl_scal, ndofsMa_scal, ndofs_scal
        integer(kind=8) :: idim, i, j, deca_row, deca_col, ndofs

        call cplCmpPenaFEMFEMMatScal(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa, &
                                     coef_pena, lhs_scal)
!
        ndim = FEFaceSl%ndim+1
        ndofsSl_scal = FEFaceSl%nbnodes
        ndofsMa_scal = FEFaceMa%nbnodes
        ndofs_scal = ndofsSl_scal+ndofsMa_scal
        ndofs = ndim*ndofs_scal
!
        lhs = 0.d0
!
        do idim = 1, ndim
            deca_col = idim
            do j = 1, ndofs_scal
                deca_row = idim
                do i = 1, j
                    lhs(deca_row, deca_col) = lhs_scal(i, j)
                    deca_row = deca_row+ndim
                end do
                deca_col = deca_col+ndim
            end do
        end do
!
! ----- Copy the lower part
!
        call hhoCopySymPartMat('U', lhs(1:ndofs, 1:ndofs))
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplCmpLagrFEMFEMMatScal(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa, FEFaceLagSl
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        real(kind=8), intent(out) :: lhs(MSIZE_LAGR_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    FEM/FEM - Compute matrix: (mu^s, u^s-u^m)
!              unknowns ((u^s, lag^s), u^m)
!===================================================================================================
!
        type(FE_Basis) :: FEBasisSl, FEBasisMa, FEBasisLagSl
        real(kind=8), dimension(MAX_BS) :: BSEvalSl, BSEvalMa, BSEvalLagSl
        real(kind=8), dimension(MAX_BS+MAX_BS) :: BS
        real(kind=8) :: coeff
        integer(kind=8) :: nbDoFsFEMa, nbDoFsFESl, nbDoFs, ipg, nbDoFsFELagSl
        blas_int, parameter :: b_one = to_blas_int(1)
        blas_int :: b_lda, b_n, b_m
!
! -- init face basis
        call FEBasisSl%initFace(FEFaceSl)
        call FEBasisMa%initFace(FEFaceMa)
        call FEBasisLagSl%initFace(FEFaceLagSl)
!
! -- number of dofs
        nbDoFsFELagSl = FEBasisLagSl%size
        nbDoFsFESl = FEBasisSl%size
        nbDoFsFEMa = FEBasisMa%size
        nbDoFs = nbDoFsFESl+nbDoFsFEMa
!
        b_n = to_blas_int(nbDoFs)
        b_m = to_blas_int(nbDoFsFELagSl)
        b_lda = to_blas_int(MSIZE_LAGR_FEM)
!
        lhs = 0.d0
!
! -- Loop on quadrature point
        do ipg = 1, FEQuadSl%nbQuadPoints
! ----- Eval face basis function at the quadrature point
! -------- Compute uh1
            BSEvalSl = FEBasisSl%func(FEQuadSl%points_param(1:3, ipg))
            BSEvalMa = FEBasisMa%func(FEQuadMa%points_param(1:3, ipg))
            BSEvalLagSl = FEBasisLagSl%func(FEQuadSl%points_param(1:3, ipg))
!
! -------- Assemble [u^s, -u^m]
            BS(1:nbDoFsFESl) = BSEvalSl(1:nbDoFsFESl)
            BS(nbDoFsFESl+1:nbDoFs) = -BSEvalMa(1:nbDoFsFEMa)
!
! --------  Eval

            coeff = FEQuadSl%weights(ipg)
            call dger(b_m, b_n, coeff, BSEvalLagSl, b_one, BS, b_one, lhs, b_lda)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine cplCmpLagrFEMFEMMatVec(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: FEFaceSl, FEFaceMa, FEFaceLagSl
        type(FE_Quadrature), intent(in) :: FEQuadSl, FEQuadMa
        real(kind=8), intent(out) :: lhs(MSIZE_LAGR_FEM, MSIZE_CPL_PENA_FEM)
!
!===================================================================================================
!    FEM/FEM - Compute matrix: (mu^s, u^s-u^m)
!              unknowns ((u^s, lag^s), u^m)
!===================================================================================================
!
        real(kind=8) :: lhs_scal(MSIZE_LAGR_FEM, MSIZE_CPL_PENA_FEM)
        integer(kind=8) :: ndim, ndofsSl_scal, ndofsMa_scal, ndofs_scal, ndofsLagSl_scal
        integer(kind=8) :: idim, i, j, deca_row, deca_col, ndofs

        call cplCmpLagrFEMFEMMatScal(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, &
                                     lhs_scal)
!
        ndim = FEFaceSl%ndim+1
        ndofsLagSl_scal = FEFaceLagSl%nbnodes
        ndofsSl_scal = FEFaceSl%nbnodes
        ndofsMa_scal = FEFaceMa%nbnodes
        ndofs_scal = ndofsSl_scal+ndofsMa_scal
        ndofs = ndim*(ndofs_scal+ndofsLagSl_scal)
!
        lhs = 0.d0
!
        do idim = 1, ndim
            deca_col = idim
            do j = 1, ndofs_scal
                deca_row = idim
                do i = 1, ndofsLagSl_scal
                    lhs(deca_row, deca_col) = lhs_scal(i, j)
                    deca_row = deca_row+ndim
                end do
                deca_col = deca_col+ndim
            end do
        end do
!
    end subroutine
!===================================================================================================
!
!===================================================================================================
!
end module
