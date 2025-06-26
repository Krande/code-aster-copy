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

module HHO_gradrec_module
!
    use HHO_type
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_basis_module
    use HHO_utils_module
    use HHO_massmat_module
    use HHO_stiffmat_module
    use HHO_geometry_module
    use HHO_algebra_module
    use HHO_matrix_module
!
    implicit none
!
    private
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
#include "blas/dger.h"
#include "blas/dposv.h"
#include "blas/dsysv.h"
!
!---------------------------------------------------------------------------------------------------
!  HHO - gradient reconstruction
!
!  This module contains all the routines to compute the gradient reconstruction for HHO methods
!
!---------------------------------------------------------------------------------------------------
!
    public :: hhoGradRecVec, hhoGradRecMat, hhoGradRecFullVec, hhoGradRecFullMat
    public :: hhoGradRecSymFullMat, hhoGradRecSymMat, hhoGradRecFullMatFromVec
!    private ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradRecVec(hhoCell, hhoData, gradrec, lhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: gradrec
        type(HHO_matrix), optional, intent(out) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the gradient reconstruction of a scalar function in Grad(P^k+1_d(T;R^d))
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   Out gradrec     : matrix of the gradient reconstruction
!   Out, option lhs : matrix (grad u, grad v) (lhs member for laplacian problem)
!
! --------------------------------------------------------------------------------------------------
! ----- Local variables
        type(HHO_Face) :: hhoFace
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_basis_face) :: hhoBasisFace
        type(HHO_quadrature) :: hhoQuad
        type(HHO_matrix) :: BG
        real(kind=8), dimension(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL) :: stiffMat, MG
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: CGradN
        real(kind=8), dimension(3, MSIZE_CELL_SCAL) :: BSCGradEval
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL), BSFEval(MSIZE_FACE_SCAL), normal(3)
        integer(kind=8) :: ipg, dimStiffMat, ifromMG, itoMG, ifromBG, itoBG, dimMG
        integer(kind=8) :: cbs, fbs, total_dofs, iface, fromFace, toFace, cell_offset
        blas_int :: b_n, b_nhrs, b_lda, b_ldb, info, b_m
        blas_int, parameter :: b_one = to_blas_int(1)
        real(kind=8) :: start, end
!
        DEBUG_TIMER(start)
!
! -- init cell basis
        call hhoBasisCell%initialize(hhoCell)
!
! -- number of dofs
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        cell_offset = total_dofs-cbs+1
!
! -- compute stiffness matrix
        dimStiffMat = hhoBasisCell%BSSize(0, hhoData%face_degree()+1)
        call hhoStiffMatCellScal(hhoCell, 0, hhoData%face_degree()+1, stiffMat)
!
! -- extract MG:  take basis functions derivatives from degree 1 to face_degree + 1
        dimMG = hhoBasisCell%BSSize(1, hhoData%face_degree()+1)
        call hhoBasisCell%BSRange(1, hhoData%face_degree()+1, ifromMG, itoMG)
!
        MG = 0.d0
        MG(1:dimMG, 1:dimMG) = stiffMat(ifromMG:itoMG, ifromMG:itoMG)
!
! -- RHS : volumetric part
        call hhoBasisCell%BSRange(0, hhoData%cell_degree(), ifromBG, itoBG)
!
        call gradrec%initialize(dimMG, total_dofs, 0.d0)
        gradrec%m(1:dimMG, cell_offset:total_dofs) = stiffMat(ifromMG:itoMG, ifromBG:itoBG)
!
        toFace = 0
! -- Loop on the faces
        do iface = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iface)
            fromFace = toFace+1
            toFace = fromFace+fbs-1
!
            call hhoBasisFace%initialize(hhoFace)
! ----- get quadrature
            call hhoQuad%GetQuadFace(hhoface, &
                                     hhoData%face_degree()+ &
                                     max(hhoData%face_degree(), hhoData%cell_degree())+1, &
                                     param=ASTER_TRUE)
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval cell basis function at the quadrature point
                call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%cell_degree(), &
                                         BSCEval)
!
! --------- Eval face basis function at the quadrature point
                call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%face_degree(), &
                                         BSFEval)
!
! --------- Eval derivative of cell basis function at the quadrature point
                call hhoBasisCell%BSEvalGrad(hhoQuad%points(1:3, ipg), 1, &
                                             hhoData%face_degree()+1, BSCGradEval)
!
! --------  Compute grad *normal
                CGradN = 0.d0
                normal = hhoNormalFaceQP(hhoFace, hhoQuad%points_param(1:2, ipg))
                b_lda = to_blas_int(3)
                b_m = to_blas_int(hhoCell%ndim)
                b_n = to_blas_int(dimMG)
                call dgemv('T', b_m, b_n, hhoQuad%weights(ipg), BSCGradEval, &
                           b_lda, normal, b_one, 0.d0, CGradN, b_one)
!
! --------  Compute (vF, grad *normal)
                b_lda = to_blas_int(dimMG)
                b_m = to_blas_int(dimMG)
                b_n = to_blas_int(fbs)
                call dger(b_m, b_n, 1.d0, CGradN, b_one, &
                          BSFEval, b_one, gradrec%m(1:dimMG, fromFace:toFace), b_lda)
!
! --------  Compute -(vT, grad *normal)
                b_lda = to_blas_int(dimMG)
                b_m = to_blas_int(dimMG)
                b_n = to_blas_int(cbs)
                call dger(b_m, b_n, -1.d0, CGradN, b_one, &
                          BSCEval, b_one, gradrec%m(1:dimMG, cell_offset:total_dofs), b_lda)
!
            end do
!
        end do
!
! - Solve the system gradrec =(MG)^-1 * BG
        if (present(lhs)) then
            call BG%initialize(dimMG, total_dofs, 0.0)
            call BG%copy(gradrec)
        end if
!
! - Verif strange bug if info neq 0 in entry
        info = 0
        b_n = to_blas_int(dimMG)
        b_nhrs = to_blas_int(total_dofs)
        b_lda = to_blas_int(MSIZE_CELL_SCAL)
        b_ldb = to_blas_int(gradrec%max_nrows)
        call dposv('U', b_n, b_nhrs, MG, b_lda, &
                   gradrec%m, b_ldb, info)
!
! - Sucess ?
        if (info .ne. 0) then
            call utmess('F', 'HHO1_4')
        end if
!
        if (present(lhs)) then
            call lhs%initialize(total_dofs, total_dofs, 0.0)
!
! ----- Compute lhs =BG**T * gradrec
            call hho_dgemm_TN(-1.d0, BG, gradrec, 0.d0, lhs)
        end if
!
        call BG%free()
!
        DEBUG_TIMER(end)
        DEBUG_TIME("Compute hhoGradRecVec", end-start)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradRecMat(hhoCell, hhoData, gradrec, gradrec_scal, lhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: gradrec
        type(HHO_matrix), intent(out) :: gradrec_scal
        type(HHO_matrix), optional, intent(out) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the gradient reconstruction of a vectorial function
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   Out gradrec     : matrix of the gradient reconstruction
!   Out, option lhs : matrix (grad u, grad v) (lhs member for laplacian problem)
!   Out gradrec_scal: matrix of the gradient reconstruction for scalar problem
!
! --------------------------------------------------------------------------------------------------
! ----- Local variables
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_matrix) :: lhs_scal
        integer(kind=8) :: gradrec_scal_row, cbs_comp, fbs_comp, faces_dofs, cbs, fbs, gbs, gbs_sym
        integer(kind=8) :: idim, ibeginGrad, iendGrad, jbeginCell, jendCell, jbeginFace, jendFace
        integer(kind=8) :: total_dofs, iFace, jbeginVec, jendVec
        real(kind=8) :: start, end
!
        DEBUG_TIMER(start)
!
! -- number of dofs
        call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs, gbs_sym)
        faces_dofs = total_dofs-cbs
!
! -- init cell basis
        call hhoBasisCell%initialize(hhoCell)
! -- dimension of stiffnes matrix for scalar problem
        gradrec_scal_row = hhoBasisCell%BSSize(1, hhoData%face_degree()+1)
!
! -- computation of the gradient reconstruction of a scalar function
        if (present(lhs)) then
            call hhoGradRecVec(hhoCell, hhoData, gradrec_scal, lhs_scal)
        else
            call hhoGradRecVec(hhoCell, hhoData, gradrec_scal)
        end if
!
! -- copy the vectorial gradient in the matrix gradient
        cbs_comp = cbs/hhoCell%ndim
        fbs_comp = fbs/hhoCell%ndim
!
        call gradrec%initialize(gbs, total_dofs, 0.d0)
! the componants are ordered direction by direction
        do idim = 1, hhoCell%ndim
! ----- copy volumetric part
            ibeginGrad = (idim-1)*gradrec_scal_row+1
            iendGrad = ibeginGrad+gradrec_scal_row-1
            jbeginCell = (idim-1)*cbs_comp+1
            jendCell = jbeginCell+cbs_comp-1
!
            gradrec%m(ibeginGrad:iendGrad, jbeginCell:jendCell) &
                = gradrec_scal%m(1:gradrec_scal_row, 1:cbs_comp)
!
! ----- copy faces part
            do iFace = 1, hhoCell%nbfaces
                jbeginFace = cbs+(iFace-1)*fbs+(idim-1)*fbs_comp+1
                jendFace = jbeginFace+fbs_comp-1
                jbeginVec = cbs_comp+(iFace-1)*fbs_comp+1
                jendVec = jbeginVec+fbs_comp-1
!
                gradrec%m(ibeginGrad:iendGrad, jbeginFace:jendFace) &
                    = gradrec_scal%m(1:gradrec_scal_row, jbeginVec:jendVec)
            end do
!
        end do
!
        if (present(lhs)) then
            call MatScal2Vec(hhoCell, hhoData, lhs_scal, lhs)
            call lhs_scal%free()
        end if
!
        DEBUG_TIMER(end)
        DEBUG_TIME("Compute hhoGradRecMat", end-start)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradRecFullVec(hhoCell, hhoData, gradrec, lhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: gradrec
        type(HHO_matrix), optional, intent(out) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the ful gradient reconstruction of a scalar function in P^k_d(T;R^d)
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   Out gradrec     : matrix of the gradient reconstruction
!   Out, option lhs : matrix (grad u, grad v) (lhs member for laplacian problem)
!
! --------------------------------------------------------------------------------------------------
! ----- Local variables
        type(HHO_Face) :: hhoFace
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_basis_face) :: hhoBasisFace
        type(HHO_quadrature) :: hhoQuad, hhoQuadCell
        type(HHO_massmat_cell) :: massMat
        type(HHO_matrix) :: BG, SOL
        real(kind=8), dimension(3, MSIZE_CELL_SCAL) :: BSCGradEval
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: BSCEval, BSGEval, VecGrad
        real(kind=8), dimension(MSIZE_FACE_SCAL) :: BSFEval
        real(kind=8) :: normal(3)
        integer(kind=8) :: cbs, fbs, total_dofs, gbs, dimMassMat
        integer(kind=8) :: ipg, ibeginBG, iendBG, ibeginSOL, iendSOL, idim
        blas_int :: b_n, b_nhrs, b_lda, b_ldb, info, b_m
        integer(kind=8) :: iface, fromFace, toFace, cell_offset
        real(kind=8) :: start, end
        blas_int, parameter :: b_one = 1
!
        DEBUG_TIMER(start)
!
! -- init cell basis
        call hhoBasisCell%initialize(hhoCell)
!
! -- number of dofs
        call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)
        cell_offset = total_dofs-cbs+1
!
! -- compute mass matrix of P^k_d(T;R)
        call massMat%compute(hhoCell, 0, hhoData%grad_degree())
        dimMassMat = massMat%nrows
!
        toFace = 0
        call gradrec%initialize(gbs, total_dofs, 0.d0)
! -- Loop on the faces
        do iface = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iface)
            fromFace = toFace+1
            toFace = fromFace+fbs-1
!
            call hhoBasisFace%initialize(hhoFace)
! ----- get quadrature
            call hhoQuad%GetQuadFace(hhoface, &
                                     hhoData%grad_degree()+max(hhoData%face_degree(), hhoData%ce&
                                     &ll_degree())+1, &
                                     param=ASTER_TRUE)
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval cell basis function at the quadrature point
                call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%cell_degree(), &
                                         BSCEval)
!
! --------- Eval face basis function at the quadrature point
                call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%face_degree(), &
                                         BSFEval)
!
! --------- Eval grad cell basis function at the quadrature point
                call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%grad_degree(), &
                                         BSGEval)
!
                normal = hhoNormalFaceQP(hhoFace, hhoQuad%points_param(1:2, ipg))
!
                do idim = 1, hhoCell%ndim
                    ibeginBG = (idim-1)*dimMassMat+1
                    iendBG = ibeginBG+dimMassMat-1
!
! ------------  Compute (vF, tau *normal)
                    b_lda = to_blas_int(dimMassMat)
                    b_m = to_blas_int(dimMassMat)
                    b_n = to_blas_int(fbs)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(idim), BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, fromFace:toFace), b_lda)
!
! ------------  Compute -(vT, tau *normal)
                    b_lda = to_blas_int(dimMassMat)
                    b_m = to_blas_int(dimMassMat)
                    b_n = to_blas_int(cbs)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(idim), BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, cell_offset:total_dofs), b_lda)
!
                end do
            end do
!
        end do
!
! -- RHS : volumetric part
! -- get quadrature
        call hhoQuadCell%GetQuadCell(hhoCell, hhoData%grad_degree()+max(hhoData%cell_degree()-1, 0))
!
! - Loop on quadrature point
        do ipg = 1, hhoQuadCell%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            call hhoBasisCell%BSEval(hhoQuadCell%points(1:3, ipg), 0, hhoData%grad_degree(), &
                                     BSGEval)
!
! ----- Eval derivative of cell basis function at the quadrature point
            call hhoBasisCell%BSEvalGrad(hhoQuadCell%points(1:3, ipg), 0, hhoData%cell_degree(), &
                                         BSCGradEval)
!
! ------ Loop on the dimension
            do idim = 1, hhoCell%ndim
                ibeginBG = (idim-1)*dimMassMat+1
                iendBG = ibeginBG+dimMassMat-1
!
                VecGrad(1:cbs) = reshape(BSCGradEval(idim, 1:cbs), (/cbs/))
! --------  Compute (grad(vT), tau)_T
                b_lda = to_blas_int(dimMassMat)
                b_m = to_blas_int(dimMassMat)
                b_n = to_blas_int(cbs)
                call dger(b_m, b_n, hhoQuadCell%weights(ipg), BSGEval, b_one, &
                          VecGrad, b_one, &
                          gradrec%m(ibeginBG:iendBG, cell_offset:total_dofs), b_lda)
            end do
        end do
!
        if (present(lhs)) then
            call BG%initialize(gbs, total_dofs, 0.0)
            call BG%copy(gradrec)
        end if
!
        if (.not. massMat%isIdentity) then
! - Solve the system gradrec =(MG)^-1 * BG
            call SOL%initialize(cbs, 3*total_dofs, 0.d0)
            do idim = 1, hhoCell%ndim
                ibeginBG = (idim-1)*dimMassMat+1
                iendBG = ibeginBG+dimMassMat-1
                ibeginSOL = (idim-1)*total_dofs+1
                iendSOL = ibeginSOL+total_dofs-1
                SOL%m(1:dimMassMat, ibeginSOL:iendSOL) = gradrec%m(ibeginBG:iendBG, 1:total_dofs)
            end do
!
! - Verif strange bug if info neq 0 in entry
            info = 0
            b_n = to_blas_int(dimMassMat)
            b_nhrs = to_blas_int(hhoCell%ndim*total_dofs)
            b_lda = to_blas_int(massMat%max_nrows)
            b_ldb = to_blas_int(SOL%max_nrows)
            call dposv('U', b_n, b_nhrs, massMat%m, b_lda, SOL%m, b_ldb, info)
!
! - Sucess ?
            if (info .ne. 0) then
                call utmess('F', 'HHO1_4')
            end if
!
! -- decompress solution
            do idim = 1, hhoCell%ndim
                ibeginBG = (idim-1)*dimMassMat+1
                iendBG = ibeginBG+dimMassMat-1
                ibeginSOL = (idim-1)*total_dofs+1
                iendSOL = ibeginSOL+total_dofs-1
                gradrec%m(ibeginBG:iendBG, 1:total_dofs) = SOL%m(1:dimMassMat, ibeginSOL:iendSOL)
            end do
            call SOL%free()
        end if
!
        if (present(lhs)) then
            call lhs%initialize(total_dofs, total_dofs, 0.d0)
!
! ----- Compute lhs =BG**T * gradrec
            call hho_dgemm_TN(1.d0, BG, gradrec, 0.0, lhs)
        end if
!
        call BG%free()
!
        DEBUG_TIMER(end)
        DEBUG_TIME("Compute hhoGradRecFullVec", end-start)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradRecFullMatFromVec(hhoCell, hhoData, gradrecvec, gradrec, lhsvec, lhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(in) :: gradrecvec
        type(HHO_matrix), intent(out) :: gradrec
        type(HHO_matrix), optional, intent(in) :: lhsvec
        type(HHO_matrix), optional, intent(out) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the ful gradient reconstruction of a vectorial function in P^k_d(T;R^(dxd))
!   In hhoCell      : the current HHO Cell
!   In hhoDta       : information on HHO methods
!   In gradrecvec   : matrix of the gradient reconstruction (for vector)
!   Out gradrec     : matrix of the gradient reconstruction (for matrix)
!   Out, option lhsvec : matrix (grad u, grad v) (lhs member for laplacian problem) (for vector)
!   Out, option lhs : matrix (grad u, grad v) (lhs member for laplacian problem) (for matrix)
!
! --------------------------------------------------------------------------------------------------
! ----- Local variables
        integer(kind=8) :: cbs_comp, fbs_comp, faces_dofs, cbs, fbs
        integer(kind=8) :: idim, ibeginGrad, iendGrad, jbeginCell, jendCell, jbeginFace, jendFace
        integer(kind=8) :: total_dofs, gbs, gbs_comp, gbs_sym, iFace, jbeginVec, jendVec
        integer(kind=8) :: faces_dofs_comp, total_dofs_comp
!
! -- number of dofs
        call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                           gbs, gbs_sym)
        faces_dofs = total_dofs-cbs
!
! -- copy the vectorial gradient in the matrix gradient
        cbs_comp = cbs/hhoCell%ndim
        fbs_comp = fbs/hhoCell%ndim
        gbs_comp = gbs/hhoCell%ndim
        faces_dofs_comp = faces_dofs/hhoCell%ndim
        total_dofs_comp = total_dofs/hhoCell%ndim
!
! -- BE CAREFULL : the componant of the gradient are stored in a row format
!
        call gradrec%initialize(gbs, total_dofs, 0.0)
        do idim = 1, hhoCell%ndim
! ----- copy volumetric part
            ibeginGrad = (idim-1)*gbs_comp+1
            iendGrad = ibeginGrad+gbs_comp-1
            jbeginCell = faces_dofs+(idim-1)*cbs_comp+1
            jendCell = jbeginCell+cbs_comp-1
!
            gradrec%m(ibeginGrad:iendGrad, jbeginCell:jendCell) &
                = gradrecvec%m(1:gbs_comp, faces_dofs_comp+1:total_dofs_comp)
!
! ----- copy faces part
            do iFace = 1, hhoCell%nbfaces
                jbeginFace = (iFace-1)*fbs+(idim-1)*fbs_comp+1
                jendFace = jbeginFace+fbs_comp-1
                jbeginVec = (iFace-1)*fbs_comp+1
                jendVec = jbeginVec+fbs_comp-1
!
                gradrec%m(ibeginGrad:iendGrad, jbeginFace:jendFace) &
                    = gradrecvec%m(1:gbs_comp, jbeginVec:jendVec)
            end do
!
        end do
!
        if (present(lhs)) then
            ASSERT(present(lhsvec))
            call lhs%initialize(total_dofs, total_dofs, 0.0)
            call MatScal2Vec(hhoCell, hhoData, lhsvec, lhs)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradRecFullMat(hhoCell, hhoData, gradrec, lhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: gradrec
        type(HHO_matrix), optional, intent(out) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the ful gradient reconstruction of a vectorial function in P^k_d(T;R^(dxd))
!   In hhoCell      : the current HHO Cell
!   In hhoDta       : information on HHO methods
!   Out gradrec     : matrix of the gradient reconstruction
!   Out, option lhs : matrix (grad u, grad v) (lhs member for laplacian problem)
!
! --------------------------------------------------------------------------------------------------
! ----- Local variables
        type(HHO_matrix) :: lhs_scal, gradrec_scal
        real(kind=8) :: start, end
!
        DEBUG_TIMER(start)
!
! -- computation of the gradient reconstruction of a scalar function
        if (present(lhs)) then
            call hhoGradRecFullVec(hhoCell, hhoData, gradrec_scal, lhs_scal)
            call hhoGradRecFullMatFromVec(hhoCell, hhoData, gradrec_scal, gradrec, &
                                          lhsvec=lhs_scal, lhs=lhs)
            call gradrec_scal%free()
            call lhs_scal%free()
        else
            call hhoGradRecFullVec(hhoCell, hhoData, gradrec_scal)
            call hhoGradRecFullMatFromVec(hhoCell, hhoData, gradrec_scal, gradrec)
            call gradrec_scal%free()
        end if
!
        DEBUG_TIMER(end)
        DEBUG_TIME("Compute hhoGradRecFullMat", end-start)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradRecSymFullMat(hhoCell, hhoData, gradrec, lhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: gradrec
        type(HHO_matrix), optional, intent(out) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Compute the full symmetric gradient reconstruction of a vect function in P^k_d(T;R^(dxd)_sym)
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   Out gradrec     : matrix of the symmetric gradient reconstruction
!   Out, option lhs : matrix (grad_s u, grad_s v) (lhs member for the symmetric laplacian problem)
!
! --------------------------------------------------------------------------------------------------
! ----- Local variables
        type(HHO_Face) :: hhoFace
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_basis_face) :: hhoBasisFace
        type(HHO_quadrature) :: hhoQuad, hhoQuadCell
        type(HHO_massmat_cell) :: massMat
        type(HHO_matrix) :: BG, SOL
        real(kind=8), dimension(6, MSIZE_CELL_VEC) :: BVCSGradEval
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL), BSFEval(MSIZE_FACE_SCAL)
        real(kind=8) :: BSGEval(MSIZE_CELL_SCAL)
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
        real(kind=8) :: coeff, normal(3)
     integer(kind=8) :: cbs, fbs, total_dofs, gbs, dimMassMat, nbdimMat, cbs_comp, fbs_comp, gbs_sym
        integer(kind=8) :: ipg, ibeginBG, iendBG, ibeginSOL, iendSOL, idim, j, iface
        integer(kind=8) :: jbegCell, jendCell, jbegFace, jendFace, faces_dofs
        blas_int :: b_n, b_nhrs, b_lda, b_ldb, info, b_m
        real(kind=8) :: start, end
        blas_int, parameter :: b_one = 1
!
        DEBUG_TIMER(start)
!
! -- init cell basis
        call hhoBasisCell%initialize(hhoCell)
!
! -- number of independant direction of the gradient
        if (hhoCell%ndim == 2) then
            nbdimMat = 3
        else if (hhoCell%ndim == 3) then
            nbdimMat = 6
        else
            ASSERT(ASTER_FALSE)
        end if
!
! -- number of dofs
        call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs, gbs_sym)
        faces_dofs = total_dofs-cbs
!
        cbs_comp = cbs/hhoCell%ndim
        fbs_comp = fbs/hhoCell%ndim
!
! -- compute mass matrix of P^k_d(T;R)
        call massMat%compute(hhoCell, 0, hhoData%grad_degree())
        dimMassMat = massMat%nrows
!
        call gradrec%initialize(gbs_sym, total_dofs, 0.d0)
!
! -- RHS : volumetric part
! -- get quadrature
        call hhoQuadCell%GetQuadCell(hhoCell, hhoData%grad_degree()+(hhoData%cell_degree()-1))
!
! -- Loop on quadrature point
        do ipg = 1, hhoQuadCell%nbQuadPoints
! ----- Eval cell basis function at the quadrature point
            call hhoBasisCell%BSEval(hhoQuadCell%points(1:3, ipg), 0, hhoData%grad_degree(), &
                                     BSGEval)
!
! ----- Eval symmetric derivative of cell basis function at the quadrature point
            call hhoBasisCell%BVEvalSymGrad(hhoQuadCell%points(1:3, ipg), 0, &
                                            hhoData%cell_degree(), BVCSGradEval)
!
! ------ Loop on diagonal terms
            do idim = 1, hhoCell%ndim
                ibeginBG = (idim-1)*dimMassMat+1
                iendBG = ibeginBG+dimMassMat-1
!
                do j = 1, cbs
                    coeff = hhoQuadCell%weights(ipg)*BVCSGradEval(idim, j)
                    b_n = to_blas_int(dimMassMat)
                    call daxpy(b_n, coeff, BSGEval, b_one, &
                               gradrec%m(ibeginBG:iendBG, faces_dofs+j), b_one)
                end do
            end do
!
! ------ Loop on extra-diagonal terms
            if (hhoCell%ndim == 3) then
                do idim = 1, hhoCell%ndim
                    ibeginBG = (3+idim-1)*dimMassMat+1
                    iendBG = ibeginBG+dimMassMat-1
!
                    do j = 1, cbs
                        coeff = hhoQuadCell%weights(ipg)*BVCSGradEval(3+idim, j)
                        b_n = to_blas_int(dimMassMat)
                        call daxpy(b_n, coeff, BSGEval, b_one, &
                                   gradrec%m(ibeginBG:iendBG, faces_dofs+j), b_one)
                    end do
                end do
            else if (hhoCell%ndim == 2) then
                ibeginBG = (3-1)*dimMassMat+1
                iendBG = ibeginBG+dimMassMat-1
!
                do j = 1, cbs
                    coeff = hhoQuadCell%weights(ipg)*BVCSGradEval(4, j)
                    b_n = to_blas_int(dimMassMat)
                    call daxpy(b_n, coeff, BSGEval, b_one, &
                               gradrec%m(ibeginBG:iendBG, faces_dofs+j), b_one)
                end do
            else
                ASSERT(ASTER_FALSE)
            end if
!
        end do
!
! -- Loop on the faces
        do iface = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iface)
!
            call hhoBasisFace%initialize(hhoFace)
! ----- get quadrature
            call hhoQuad%GetQuadFace(hhoface, &
                                     hhoData%grad_degree()+max(hhoData%face_degree(), hhoData%ce&
                                     &ll_degree())+1, &
                                     param=ASTER_TRUE)
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval cell basis function at the quadrature point
                call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%cell_degree(), &
                                         BSCEval)
!
! --------- Eval face basis function at the quadrature point
                call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%face_degree(), &
                                         BSFEval)
!
! --------- Eval grad cell basis function at the quadrature point
                call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%grad_degree(), &
                                         BSGEval)
!
                normal = hhoNormalFaceQP(hhoFace, hhoQuad%points_param(1:2, ipg))
!
! ----- copy by dimension
! ----- diagonal composants
                do idim = 1, hhoCell%ndim
                    ibeginBG = (idim-1)*dimMassMat+1
                    iendBG = ibeginBG+dimMassMat-1
                    jbegCell = faces_dofs+(idim-1)*cbs_comp+1
                    jendCell = jbegCell+cbs_comp-1
                    jbegFace = (iface-1)*fbs+(idim-1)*fbs_comp+1
                    jendFace = jbegFace+fbs_comp-1
!
! ------------  Compute -(vT, tau *normal)
                    b_lda = to_blas_int(dimMassMat)
                    b_m = to_blas_int(dimMassMat)
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(idim), BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
! ------------  Compute (vF, tau *normal)
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(idim), BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
                end do
!
                if (hhoCell%ndim == 2) then
! ------ extra diagonal composants term 12
                    ibeginBG = 2*dimMassMat+1
                    iendBG = ibeginBG+dimMassMat-1
!
! ------------  Compute -(vT, tau *normal)
                    jbegCell = faces_dofs+1
                    jendCell = jbegCell+cbs_comp-1
                    b_lda = to_blas_int(dimMassMat)
                    b_m = to_blas_int(dimMassMat)
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(2)/rac2, BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
                    jbegCell = jbegCell+cbs_comp
                    jendCell = jbegCell+cbs_comp-1
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(1)/rac2, BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
! ------------  Compute (vF, tau *normal)
                    jbegFace = (iface-1)*fbs+1
                    jendFace = jbegFace+fbs_comp-1
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(2)/rac2, BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
!
                    jbegFace = jbegFace+fbs_comp
                    jendFace = jbegFace+fbs_comp-1
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(1)/rac2, BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
!
                else if (hhoCell%ndim == 3) then
! ------ hors diagonal composants term 12
                    ibeginBG = 3*dimMassMat+1
                    iendBG = ibeginBG+dimMassMat-1
!
! ------------  Compute -(vT, tau *normal)
                    jbegCell = faces_dofs+1
                    jendCell = jbegCell+cbs_comp-1
                    b_lda = to_blas_int(dimMassMat)
                    b_m = to_blas_int(dimMassMat)
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(2)/rac2, BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
                    jbegCell = jbegCell+cbs_comp
                    jendCell = jbegCell+cbs_comp-1
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(1)/rac2, BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
! ------------  Compute (vF, tau *normal)
                    jbegFace = (iface-1)*fbs+1
                    jendFace = jbegFace+fbs_comp-1
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(2)/rac2, BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
!
                    jbegFace = jbegFace+fbs_comp
                    jendFace = jbegFace+fbs_comp-1
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(1)/rac2, BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
! ------ extra diagonal composants term 13
                    ibeginBG = 4*dimMassMat+1
                    iendBG = ibeginBG+dimMassMat-1
!
! ------------  Compute -(vT, tau *normal)
                    jbegCell = faces_dofs+1
                    jendCell = jbegCell+cbs_comp-1
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(3)/rac2, BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
                    jbegCell = faces_dofs+2*cbs_comp+1
                    jendCell = jbegCell+cbs_comp-1
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(1)/rac2, BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
! ------------  Compute (vF, tau *normal)
                    jbegFace = (iface-1)*fbs+1
                    jendFace = jbegFace+fbs_comp-1
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(3)/rac2, BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
!
                    jbegFace = (iface-1)*fbs+1+2*fbs_comp
                    jendFace = jbegFace+fbs_comp-1
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(1)/rac2, BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
! ------ extra diagonal composants term 23
                    ibeginBG = 5*dimMassMat+1
                    iendBG = ibeginBG+dimMassMat-1
!
! ------------  Compute -(vT, tau *normal)
                    jbegCell = faces_dofs+cbs_comp+1
                    jendCell = jbegCell+cbs_comp-1
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(3)/rac2, BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
                    jbegCell = faces_dofs+2*cbs_comp+1
                    jendCell = jbegCell+cbs_comp-1
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -hhoQuad%weights(ipg)*normal(2)/rac2, BSGEval, b_one, &
                              BSCEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegCell:jendCell), b_lda)
!
! ------------  Compute (vF, tau *normal)
                    jbegFace = (iface-1)*fbs+1+fbs_comp
                    jendFace = jbegFace+fbs_comp-1
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(3)/rac2, BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
!
                    jbegFace = (iface-1)*fbs+1+2*fbs_comp
                    jendFace = jbegFace+fbs_comp-1
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, hhoQuad%weights(ipg)*normal(2)/rac2, BSGEval, b_one, &
                              BSFEval, b_one, &
                              gradrec%m(ibeginBG:iendBG, jbegFace:jendFace), b_lda)
!
                else
                    ASSERT(ASTER_FALSE)
                end if
!
            end do
!
        end do
!
        if (present(lhs)) then
            call BG%initialize(gbs_sym, total_dofs)
            call BG%copy(gradrec)
        end if
!
        if (.not. massMat%isIdentity) then
! - Solve the system gradrec =(MG)^-1 * BG
            call SOL%initialize(dimMassMat, nbdimMat*total_dofs)
            do idim = 1, nbdimMat
                ibeginBG = (idim-1)*dimMassMat+1
                iendBG = ibeginBG+dimMassMat-1
                ibeginSOL = (idim-1)*total_dofs+1
                iendSOL = ibeginSOL+total_dofs-1
                SOL%m(1:dimMassMat, ibeginSOL:iendSOL) = gradrec%m(ibeginBG:iendBG, 1:total_dofs)
            end do
!
! - Verif strange bug if info neq 0 in entry
            info = 0
            b_n = to_blas_int(dimMassMat)
            b_nhrs = to_blas_int(nbdimMat*total_dofs)
            b_lda = to_blas_int(MSIZE_CELL_SCAL)
            b_ldb = to_blas_int(SOL%max_nrows)
            call dposv('U', b_n, b_nhrs, massMat%m, b_lda, &
                       SOL%m, b_ldb, info)
!
! - Sucess ?
            if (info .ne. 0) then
                call utmess('F', 'HHO1_4')
            end if
!
! -- decompress solution
            do idim = 1, nbdimMat
                ibeginBG = (idim-1)*dimMassMat+1
                iendBG = ibeginBG+dimMassMat-1
                ibeginSOL = (idim-1)*total_dofs+1
                iendSOL = ibeginSOL+total_dofs-1
                gradrec%m(ibeginBG:iendBG, 1:total_dofs) = SOL%m(1:dimMassMat, ibeginSOL:iendSOL)
            end do
            call SOL%free()
        end if
!
        if (present(lhs)) then
            call lhs%initialize(total_dofs, total_dofs, 0.0)
!
! ----- Compute lhs =BG**T * gradrec
            call hho_dgemm_TN(1.d0, BG, gradrec, 0.d0, lhs)
!
        end if
        call BG%free()
!
        DEBUG_TIMER(end)
        DEBUG_TIME("Compute hhoGradRecSymFullMat", end-start)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradRecSymMat(hhoCell, hhoData, gradrec, lhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: gradrec
        type(HHO_matrix), optional, intent(out) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the symmetric gradient reconstruction of a vectoriel function
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   Out gradrec     : matrix of the symmetric gradient reconstruction
!   Out, option lhs : matrix (grad_s u, grad_s v) (lhs member for symmetric laplacian problem)
!
! --------------------------------------------------------------------------------------------------
! ----- Local variables
        type(HHO_Face) :: hhoFace
        type(HHO_basis_cell) :: hhoBasisCell
        type(HHO_basis_face) :: hhoBasisFace
        type(HHO_quadrature) :: hhoQuad, hhoQuadCell
        type(HHO_matrix) :: BG, stiffMat, MG
        real(kind=8), dimension(MSIZE_CELL_VEC, 3) :: CGradN
        real(kind=8), dimension(6, MSIZE_CELL_VEC) :: BVCGradEval
        real(kind=8), dimension(3, MSIZE_CELL_SCAL) :: BSCGradEval
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL), BSFEval(MSIZE_FACE_SCAL)
        type(HHO_matrix) :: WORK
        blas_int, dimension(MSIZE_CELL_VEC+3) :: IPIV
        integer(kind=8) :: ipg, dimStiffMat, ifromBG, itoBG, dimMG, nblag, dimMGLag, idir, idir2
        integer(kind=8) :: cbs, fbs, total_dofs, iface, i, cbs_comp, fbs_comp, dimMG_cmp, ind_MG
        integer(kind=8) :: jbeginCell, jendCell, jbeginFace, jendFace, idim, j, dimStiffMat_cmp
        integer(kind=8) :: row_deb_MG, row_fin_MG, col_deb_MG, col_fin_MG, col_deb_BG, col_fin_BG
        integer(kind=8) :: row_deb_ST, row_fin_ST, col_deb_ST, col_fin_ST, faces_dof
        blas_int :: b_n, b_nhrs, b_lda, b_ldb, info, LWORK, b_m
        real(kind=8) :: qp_dphi_ss, normal(3)
        real(kind=8) :: start, end
        blas_int, parameter :: b_one = 1
!
        DEBUG_TIMER(start)
!
! -- init cell basis
        call hhoBasisCell%initialize(hhoCell)
!
! -- number of dofs
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        faces_dof = total_dofs-cbs
        cbs_comp = cbs/hhoCell%ndim
        fbs_comp = fbs/hhoCell%ndim
!
! -- number of lagrange to impose skew symmetric part to zero
        nblag = 0
        if (hhoCell%ndim == 2) then
            nblag = 1
        else if (hhoCell%ndim == 3) then
            nblag = 3
        else
            ASSERT(ASTER_FALSE)
        end if
!
! -- compute stiffness matrix
        dimStiffMat = hhoBasisCell%BVSize(0, hhoData%face_degree()+1)
        dimStiffMat_cmp = dimStiffMat/hhoCell%ndim
        call hhoSymStiffMatCellVec(hhoCell, 0, hhoData%face_degree()+1, stiffMat)
!
! -- extract MG:  take basis functions derivatives from degree 1 to face_degree + 1
        dimMG = hhoBasisCell%BVSize(1, hhoData%face_degree()+1)
        dimMG_cmp = dimMG/hhoCell%ndim
        dimMGLag = dimMG+nblag
!
        call MG%initialize(dimMG+3, dimMG+3, 0.d0)
        call BG%initialize(dimMGLag, total_dofs, 0.d0)
!
        do idir = 1, hhoCell%ndim
            row_deb_MG = (idir-1)*dimMG_cmp+1
            row_fin_MG = row_deb_MG+dimMG_cmp-1
            row_deb_ST = (idir-1)*dimStiffMat_cmp+2
            row_fin_ST = row_deb_ST+dimMG_cmp-1
!
            do idir2 = 1, hhoCell%ndim
                col_deb_MG = (idir2-1)*dimMG_cmp+1
                col_fin_MG = col_deb_MG+dimMG_cmp-1
                col_deb_ST = (idir2-1)*dimStiffMat_cmp+2
                col_fin_ST = col_deb_ST+dimMG_cmp-1
!
                MG%m(row_deb_MG:row_fin_MG, col_deb_MG:col_fin_MG) = stiffMat%m( &
                                                                     row_deb_ST:row_fin_ST, &
                                                                     col_deb_ST:col_fin_ST &
                                                                     )
!
                MG%m(col_deb_MG:col_fin_MG, row_deb_MG:row_fin_MG) = stiffMat%m( &
                                                                     col_deb_ST:col_fin_ST, &
                                                                     row_deb_ST:row_fin_ST &
                                                                     )
!
                col_deb_BG = faces_dof+(idir2-1)*cbs_comp+1
                col_fin_BG = col_deb_BG+cbs_comp-1
                col_deb_ST = (idir2-1)*dimStiffMat_cmp+1
                col_fin_ST = col_deb_ST+cbs_comp-1
!
                BG%m(row_deb_MG:row_fin_MG, col_deb_BG:col_fin_BG) = stiffMat%m( &
                                                                     row_deb_ST:row_fin_ST, &
                                                                     col_deb_ST:col_fin_ST &
                                                                     )
            end do
!
        end do
!
        call stiffMat%free()
!
! -- impose lagrange multipliers
! -- get quadrature
        call hhoQuadCell%GetQuadCell(hhoCell, hhoData%face_degree()+1)
!
! -- Loop on quadrature point
        do ipg = 1, hhoQuadCell%nbQuadPoints
!
! ----- Eval derivative of cell basis function at the quadrature point
            call hhoBasisCell%BSEvalGrad(hhoQuadCell%points(1:3, ipg), 1, &
                                         hhoData%face_degree()+1, BSCGradEval)
            do j = 1, dimMG_cmp
                if (hhoCell%ndim == 2) then
! ------------- lag1
! ------------- dir = 1
                    qp_dphi_ss = hhoQuadCell%weights(ipg)*BSCGradEval(2, j)
                    ind_MG = j
                    MG%m(ind_MG, dimMG+1) = MG%m(ind_MG, dimMG+1)+qp_dphi_ss
                    MG%m(dimMG+1, ind_MG) = MG%m(dimMG+1, ind_MG)+qp_dphi_ss
! ------------- dir = 2
                    qp_dphi_ss = -hhoQuadCell%weights(ipg)*BSCGradEval(1, j)
                    ind_MG = dimMG_cmp+j
                    MG%m(ind_MG, dimMG+1) = MG%m(ind_MG, dimMG+1)+qp_dphi_ss
                    MG%m(dimMG+1, ind_MG) = MG%m(dimMG+1, ind_MG)+qp_dphi_ss
                else if (hhoCell%ndim == 3) then
! ------------- lag1
! ------------- dir = 1
                    qp_dphi_ss = hhoQuadCell%weights(ipg)*BSCGradEval(2, j)
                    ind_MG = j
                    MG%m(ind_MG, dimMG+1) = MG%m(ind_MG, dimMG+1)+qp_dphi_ss
                    MG%m(dimMG+1, ind_MG) = MG%m(dimMG+1, ind_MG)+qp_dphi_ss
! ------------- dir = 2
                    qp_dphi_ss = -hhoQuadCell%weights(ipg)*BSCGradEval(1, j)
                    ind_MG = dimMG_cmp+j
                    MG%m(ind_MG, dimMG+1) = MG%m(ind_MG, dimMG+1)+qp_dphi_ss
                    MG%m(dimMG+1, ind_MG) = MG%m(dimMG+1, ind_MG)+qp_dphi_ss
! ------------- lag2
! ------------- dir = 1
                    qp_dphi_ss = hhoQuadCell%weights(ipg)*BSCGradEval(3, j)
                    ind_MG = j
                    MG%m(ind_MG, dimMG+2) = MG%m(ind_MG, dimMG+2)+qp_dphi_ss
                    MG%m(dimMG+2, ind_MG) = MG%m(dimMG+2, ind_MG)+qp_dphi_ss
! ------------- dir = 3
                    qp_dphi_ss = -hhoQuadCell%weights(ipg)*BSCGradEval(1, j)
                    ind_MG = 2*dimMG_cmp+j
                    MG%m(ind_MG, dimMG+2) = MG%m(ind_MG, dimMG+2)+qp_dphi_ss
                    MG%m(dimMG+2, ind_MG) = MG%m(dimMG+2, ind_MG)+qp_dphi_ss
! ------------- lag3
! ------------- dir = 2
                    qp_dphi_ss = hhoQuadCell%weights(ipg)*BSCGradEval(3, j)
                    ind_MG = dimMG_cmp+j
                    MG%m(ind_MG, dimMG+3) = MG%m(ind_MG, dimMG+3)+qp_dphi_ss
                    MG%m(dimMG+3, ind_MG) = MG%m(dimMG+3, ind_MG)+qp_dphi_ss
! ------------- dir = 3
                    qp_dphi_ss = -hhoQuadCell%weights(ipg)*BSCGradEval(2, j)
                    ind_MG = 2*dimMG_cmp+j
                    MG%m(ind_MG, dimMG+3) = MG%m(ind_MG, dimMG+3)+qp_dphi_ss
                    MG%m(dimMG+3, ind_MG) = MG%m(dimMG+3, ind_MG)+qp_dphi_ss
                else
                    ASSERT(ASTER_FALSE)
                end if
            end do
        end do
!
! -- RHS : volumetric part
        call hhoBasisCell%BVRange(0, hhoData%cell_degree(), ifromBG, itoBG)
!
! -- Loop on the faces
        do iface = 1, hhoCell%nbfaces
            hhoFace = hhoCell%faces(iface)
!
            call hhoBasisFace%initialize(hhoFace)
! ----- get quadrature
            call hhoQuad%GetQuadFace(hhoface, &
                                     hhoData%face_degree()+max(hhoData%face_degree(), hhoData%ce&
                                     &ll_degree()), &
                                     param=ASTER_TRUE)
!
! ----- Loop on quadrature point
            do ipg = 1, hhoQuad%nbQuadPoints
! --------- Eval cell basis function at the quadrature point
                call hhoBasisCell%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%cell_degree(), &
                                         BSCEval)
!
! --------- Eval face basis function at the quadrature point
                call hhoBasisFace%BSEval(hhoQuad%points(1:3, ipg), 0, hhoData%face_degree(), &
                                         BSFEval)
!
! --------- Eval symetric derivative of cell basis function at the quadrature point
                call hhoBasisCell%BVEvalSymGrad(hhoQuad%points(1:3, ipg), 1, &
                                                hhoData%face_degree()+1, BVCGradEval)
!
! --------  Compute grad_s *normal
                CGradN = 0.d0
                normal = hhoNormalFaceQP(hhoFace, hhoQuad%points_param(1:2, ipg))
                do i = 1, dimMG
                    call hhoProdSmatVec(BVCGradEval(1:6, i), normal, hhoCell%ndim, &
                                        CGradN(i, 1:3))
                    CGradN(i, 1:3) = hhoQuad%weights(ipg)*CGradN(i, 1:3)
                end do
!
                do idim = 1, hhoCell%ndim
! ------------- Compute (vF, grad_s *normal)
                    jbeginFace = (idim-1)*fbs_comp+(iface-1)*fbs+1
                    jendFace = jbeginFace+fbs_comp-1
!
                    b_lda = to_blas_int(dimMG)
                    b_m = to_blas_int(dimMG)
                    b_n = to_blas_int(fbs_comp)
                    call dger(b_m, b_n, 1.d0, CGradN(1:dimMG, idim), b_one, &
                              BSFEval, b_one, BG%m(1:dimMG, jbeginFace:jendFace), b_lda)
!
! ------------  Compute -(vT, grad_s *normal)
                    jbeginCell = faces_dof+(idim-1)*cbs_comp+1
                    jendCell = jbeginCell+cbs_comp-1
                    b_lda = to_blas_int(dimMG)
                    b_m = to_blas_int(dimMG)
                    b_n = to_blas_int(cbs_comp)
                    call dger(b_m, b_n, -1.d0, CGradN(1:dimMG, idim), b_one, &
                              BSCEval, b_one, BG%m(1:dimMG, jbeginCell:jendCell), b_lda)
                end do
!
            end do
!
        end do
!
! - Solve the system gradrec =(MG)^-1 * BG
        call gradrec%initialize(dimMGLag, total_dofs)
        call gradrec%copy(BG)
!
! - Verif strange bug if info neq 0 in entry
        call WORK%initialize(dimMGLag, total_dofs)
        LWORK = to_blas_int(dimMGLag*total_dofs)
        info = 0
        b_n = to_blas_int(dimMGLag)
        b_nhrs = to_blas_int(total_dofs)
        b_lda = to_blas_int(MG%max_nrows)
        b_ldb = to_blas_int(gradrec%max_nrows)
        call dsysv('U', b_n, b_nhrs, MG%m, b_lda, &
                   IPIV, gradrec%m, b_ldb, WORK%m, LWORK, info)
        call WORK%free()
!
! - Resize matrix to remove lagr
        gradrec%nrows = dimMG
!
! - Sucess ?
        if (info .ne. 0) then
            call utmess('F', 'HHO1_5')
        end if
!
        if (present(lhs)) then
            call lhs%initialize(total_dofs, total_dofs, 0.0)
!
! ----- Compute lhs =BG**T * gradrec
            call hho_dgemm_TN(1.d0, BG, gradrec, 0.d0, lhs)
        end if

        call BG%free()
        call MG%free()
!
        DEBUG_TIMER(end)
        DEBUG_TIME("Compute hhoGradRecSymMat", end-start)
!
    end subroutine
!
end module
