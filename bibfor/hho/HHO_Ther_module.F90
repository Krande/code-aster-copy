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

module HHO_Ther_module
!
    use FE_algebra_module
    use HHO_basis_module
    use HHO_Dirichlet_module
    use HHO_eval_module
    use HHO_gradrec_module, only: hhoGradRecFullVec, hhoGradRecVec
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_stabilization_module, only: hhoStabScal, hdgStabScal
    use HHO_type
    use HHO_utils_module
    use HHO_matrix_module
    use HHO_algebra_module
    use NonLin_Datastructure_type
!
    implicit none
!
    private
#include "asterf_debug.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/matrRotLGTher.h"
#include "asterfort/nlcomp.h"
#include "asterfort/ntfcma.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcfode.h"
#include "asterfort/rcvalb.h"
#include "asterfort/readMatrix.h"
#include "asterfort/readVector.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dgemm.h"
#include "blas/dger.h"
#include "blas/dsymv.h"
#include "blas/dsyr.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - thermics
!
! Specific routines for thermics
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: hhoLocalRigiTher, hhoCalcStabCoeffTher, hhoLocalMassTher
    public :: hhoCalcOpTher, hhoComputeRhsRigiTher, hhoComputeLhsRigiTher
    public :: hhoComputeLhsMassTher, hhoComputeRhsMassTher, hhoComputeBehaviourTher
    public :: hhoReloadPreCalcTher, hhoLocalMassTherNL
    private :: hhoComputeAgphi
    private :: LambdaMax, hhoComputeRhoCpTher, hhoComputeBetaTher
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLocalRigiTher(hhoCell, hhoData, hhoQuadCellRigi, option, gradrec, &
                                stab, fami, lhs, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellRigi
        character(len=16), intent(in) :: option
        type(HHO_matrix), intent(in) :: gradrec
        type(HHO_matrix), intent(in) :: stab
        character(len=8), intent(in) :: fami
        type(HHO_matrix), intent(out), optional :: lhs
        real(kind=8), intent(out), optional :: rhs(MSIZE_TDOFS_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the local rigidity contribution for thermics
!   RHS = -(lambda * \GkT(huT), \GkT(hvT))_T + lambda * s_T(huT, hvT)
!   LHS = (lambda * \GkT(hduT), \GkT(hvT))_T + lambda * s_T(hduT, hvT)
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   In gradrec      : local gradient reconstruction
!   In stab         : local stabilisation
!   In fami         : familly of quadrature points (of hhoQuadCellRigi)
!   In typmod       : type of modelization
!   In compor       : type of behavior
!   In option       : option of computations
!   Out lhs         : local contribution (lhs)
!   Out rhs         : local contribution (rhs)
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        character(len=32) :: phenom
        integer :: cbs, fbs, total_dofs, j, faces_dofs, gbs
        integer :: jmate, ipg, icodre(3), jtemps, ndim
        real(kind=8), dimension(MSIZE_CELL_VEC) :: bT, G_curr_coeff
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: temp_curr
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8) :: AT(MSIZE_CELL_VEC, MSIZE_CELL_VEC)
        type(HHO_matrix) :: TMP
        real(kind=8) :: module_tang(3, 3), G_curr(3), sig_curr(3)
        real(kind=8) :: coorpg(3), weight, time_curr, temp_eval_curr
        real(kind=8), pointer :: flux(:) => null()
        aster_logical :: l_rhs, l_lhs, l_nl, l_flux
        real(kind=8) :: start, end
        blas_int :: b_incx, b_incy, b_n
        blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m
!
        DEBUG_TIMER(start)
!
        l_lhs = present(lhs)
        l_rhs = present(rhs)
!
        l_nl = ASTER_FALSE
        if (option == "RIGI_THER_TANG" .or. option == "RAPH_THER") then
            l_nl = ASTER_TRUE
        end if
        l_flux = ASTER_FALSE
        if (option == "RAPH_THER") then
            l_flux = ASTER_TRUE
            call jevech("PFLUXPR", "E", vr=flux)
        end if
!
! --- Get input fields
!
        call jevech('PMATERC', 'L', jmate)
!
        call rccoma(zi(jmate), 'THER', 1, phenom, icodre(1))
!
        time_curr = 0.d0
        if (.not. l_nl) then
            call jevech('PINSTR', 'L', jtemps)
            time_curr = zr(jtemps)
        end if
!
! --- number of dofs
!
        call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)
        faces_dofs = total_dofs-cbs
!
! -- initialization
!
        if (l_lhs) then
            call lhs%initialize(total_dofs, total_dofs, 0.d0)
        end if
        if (l_rhs) rhs = 0.d0
!
        bT = 0.d0
        AT = 0.d0
        G_curr_coeff = 0.d0
!
        call hhoBasisCell%initialize(hhoCell)
!
! --- compute temp in T+
!
        temp_curr = 0.d0
        if (l_nl) then
            call readVector('PTEMPEI', total_dofs, temp_curr)
            call hhoRenumTherVecInv(hhoCell, hhoData, temp_curr)
        else if (l_rhs) then
            call readVector('PTEMPER', total_dofs, temp_curr)
            call hhoRenumTherVecInv(hhoCell, hhoData, temp_curr)
        end if
        ndim = hhoCell%ndim
!
! ----- compute G_curr = gradrec * temp_curr
!
        call hho_dgemv_N(1.d0, gradrec, temp_curr, 0.0, G_curr_coeff)
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCellRigi%nbQuadPoints
            coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
            weight = hhoQuadCellRigi%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
            call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
! --------- Eval gradient at T+
!
            G_curr = hhoEvalVecCell( &
                     hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), G_curr_coeff, gbs)
!
! --------- Eval temperature at T+
!
            temp_eval_curr = hhoEvalScalCell( &
                             hhoBasisCell, hhoData%cell_degree(), coorpg(1:3), temp_curr, cbs)
!
! ------- Compute behavior
!
            call hhoComputeBehaviourTher(phenom, fami, ipg, ndim, time_curr, &
                                         jmate, coorpg, temp_eval_curr, G_curr, sig_curr, &
                                         module_tang)
!
! ------- Compute rhs
!
            if (l_rhs) call hhoComputeRhsRigiTher(hhoCell, sig_curr, weight, BSCEval, gbs, &
                                                  bT)
!
! ------- Compute lhs
!
            if (l_lhs) call hhoComputeLhsRigiTher(hhoCell, module_tang, weight, BSCEval, gbs, &
                                                  AT)
!
! ------- Save fluxes
!
            if (l_flux) then
                flux(ndim*(ipg-1)+1:ndim*(ipg-1)+ndim) = -sig_curr(1:ndim)
            end if
!
        end do
!
! ----- compute rhs += Gradrec**T * bT
!
        if (l_rhs) then
            call hho_dgemv_T(1.0, gradrec, bT, 1.d0, rhs)
        end if
!
! ----- compute lhs += gradrec**T * AT * gradrec
! ----- step1: TMP = AT * gradrec
!
        if (l_lhs) then
            call TMP%initialize(gbs, total_dofs, 0.d0)
            b_ldc = to_blas_int(TMP%max_nrows)
            b_ldb = to_blas_int(gradrec%max_nrows)
            b_lda = to_blas_int(MSIZE_CELL_VEC)
            b_m = to_blas_int(gbs)
            b_n = to_blas_int(total_dofs)
            b_k = to_blas_int(gbs)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       1.d0, AT, b_lda, gradrec%m, b_ldb, &
                       0.d0, TMP%m, b_ldc)
!
! ----- step2: lhs += gradrec**T * TMP
!
            call hho_dgemm_TN(1.d0, gradrec, TMP, 1.d0, lhs)
            call TMP%free()
!
        end if
!
! --- add stabilization
!
        call hhoCalcStabCoeffTher(hhoData, fami, hhoQuadCellRigi%nbQuadPoints, hhoCell%ndim, &
                                  time_curr)
!
        if (l_rhs) then
            b_lda = to_blas_int(stab%max_nrows)
            b_n = to_blas_int(total_dofs)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dsymv('U', b_n, hhoData%coeff_stab(), stab%m, b_lda, &
                       temp_curr, b_incx, 1.d0, rhs, b_incy)
        end if
!
        if (l_lhs) then
            do j = 1, total_dofs
                b_n = to_blas_int(total_dofs)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, hhoData%coeff_stab(), stab%m(:, j), b_incx, lhs%m(:, j), &
                           b_incy)
            end do
        end if
!
        DEBUG_TIMER(end)
        DEBUG_TIME("Compute hhoLocalRigiTher ("//option//")", end-start)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLocalMassTher(hhoCell, hhoData, hhoQuadCellMass, fami, lhs, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellMass
        character(len=8), intent(in) :: fami
        type(HHO_matrix), intent(out), optional :: lhs
        real(kind=8), intent(out), optional :: rhs(MSIZE_TDOFS_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the local mass contribution for thermics
!   RHS = (rho_cp * uT, vT)_T
!   LHS = (rho_cp * duT, vT)_T
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   In fami         : familly of quadrature points (of hhoQuadCellRigi)
!   In typmod       : type of modelization
!   In compor       : type of behavior
!   In option       : option of computations
!   Out lhs         : local contribution (lhs)
!   Out rhs         : local contribution (rhs)
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        character(len=32) :: phenom
        integer :: cbs, fbs, total_dofs, faces_dofs, gbs
        integer :: jmate, ipg, icodre(3), jtemps
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: temp_T_curr
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8) :: coorpg(3), weight, time_curr, cp, temp_eval, beta
        character(len=8) :: poum
        aster_logical :: l_rhs, l_lhs
!
        l_lhs = present(lhs)
        l_rhs = present(rhs)
!
! --- Get input fields
!
        call jevech('PMATERC', 'L', jmate)
!
        call rccoma(zi(jmate), 'THER', 1, phenom, icodre(1))
!
        call jevech('PINSTR', 'L', jtemps)
        time_curr = zr(jtemps)
!
! --- number of dofs
!
        call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                           gbs)
        faces_dofs = total_dofs-cbs
!
! -- initialization
!
        if (l_lhs) then
            call lhs%initialize(total_dofs, total_dofs, 0.d0)
        end if
        if (l_rhs) rhs = 0.d0
!
        call hhoBasisCell%initialize(hhoCell)
!
! --- compute temp in T+
!
        temp_T_curr = 0.d0
        if (l_rhs) call readVector('PTEMPER', cbs, temp_T_curr, total_dofs-cbs)
        poum = "+"
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCellMass%nbQuadPoints
            coorpg(1:3) = hhoQuadCellMass%points(1:3, ipg)
            weight = hhoQuadCellMass%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
            call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%cell_degree(), BSCEval)
!
! --------- Eval gradient at T+
!
            temp_eval = hhoEvalScalCell( &
                        hhoBasisCell, hhoData%cell_degree(), coorpg, temp_T_curr, cbs)
!
! -------- Compute behavior
!
            call hhoComputeRhoCpTher(phenom, fami, poum, ipg, time_curr, &
                                     jmate, cp)
!
! -------- Compute rhs
!
            if (l_rhs) then
                beta = cp*temp_eval
                call hhoComputeRhsMassTher(beta, weight, BSCEval, cbs, rhs(1:MSIZE_CELL_SCAL))
            end if
!
! -------- Compute lhs
!
            if (l_lhs) then
                call hhoComputeLhsMassTher(cp, weight, BSCEval, cbs, lhs)
            end if
!
        end do
!
! ----- Copy the lower part
!
        if (l_lhs) call hhoCopySymPartMat('U', lhs%m(1:cbs, 1:cbs))
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLocalMassTherNL(hhoCell, hhoData, hhoQuadCellMass, lhs, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellMass
        type(HHO_matrix), intent(out), optional :: lhs
        real(kind=8), intent(out), optional :: rhs(MSIZE_TDOFS_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the nonlinear local mass contribution for thermics
!   RHS = (beta(T) * uT, vT)_T
!   LHS = (dbeta(T)_dT * duT, vT)_T
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellMas : quadrature rules from the rigidity family
!   In fami         : familly of quadrature points (of hhoQuadCellRigi)
!   Out lhs         : local contribution (lhs)
!   Out rhs         : local contribution (rhs)
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        character(len=32) :: phenom
        integer :: cbs, fbs, total_dofs, faces_dofs, gbs
        integer :: jmate, ipg, icodre(3), jcomp, ifon(6)
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: temp_T_curr
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8) :: coorpg(3), weight, temp_eval, beta, dbeta
        character(len=16) :: comp
        aster_logical :: l_rhs, l_lhs, aniso
!
        l_lhs = present(lhs)
        l_rhs = present(rhs)
!
! --- Get input fields
!
        call jevech('PMATERC', 'L', jmate)
        call jevech('PCOMPOR', 'L', jcomp)
!
        comp = zk16(jcomp)
        if (comp(1:5) .eq. 'THER_') then
!
            call rccoma(zi(jmate), 'THER', 1, phenom, icodre(1))
            aniso = ASTER_FALSE
            if (phenom(1:12) .eq. 'THER_NL_ORTH') then
                aniso = ASTER_TRUE
            end if
            call ntfcma(comp, zi(jmate), aniso, ifon)
            if (comp(1:9) .eq. 'THER_HYDR') then
                ASSERT(ASTER_FALSE)
            end if
        end if
!
! --- number of dofs
!
        call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                           gbs)
        faces_dofs = total_dofs-cbs
!
! -- initialization
!
        if (l_lhs) then
            call lhs%initialize(total_dofs, total_dofs, 0.0)
        end if
        if (l_rhs) rhs = 0.d0
!
        call hhoBasisCell%initialize(hhoCell)
!
! --- compute temp in T+
!
        temp_T_curr = 0.d0
        call readVector('PTEMPEI', cbs, temp_T_curr, total_dofs-cbs)
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCellMass%nbQuadPoints
            coorpg(1:3) = hhoQuadCellMass%points(1:3, ipg)
            weight = hhoQuadCellMass%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
            call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%cell_degree(), BSCEval)
!
! --------- Eval gradient at T+
!
            temp_eval = hhoEvalScalCell( &
                        hhoBasisCell, hhoData%cell_degree(), coorpg, temp_T_curr, cbs)
!
! -------- Compute behavior
!
            call hhoComputeBetaTher(comp, ifon(1), temp_eval, beta, dbeta)
!
! -------- Compute rhs
!
            if (l_rhs) then
                call hhoComputeRhsMassTher(beta, weight, BSCEval, cbs, rhs(1:MSIZE_CELL_SCAL))
            end if
!
! -------- Compute lhs
!
            if (l_lhs) then
                call hhoComputeLhsMassTher(dbeta, weight, BSCEval, cbs, lhs)
            end if
!
        end do
!
! ----- Copy the lower part
!
        if (l_lhs) call hhoCopySymPartMat('U', lhs%m(1:cbs, 1:cbs))
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoReloadPreCalcTher(hhoCell, hhoData, gradfull, stab)
!
        implicit none
!
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_matrix), intent(out) :: gradfull
        type(HHO_matrix), optional, intent(out) :: stab
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - Reload Precomputation of operators
!
! In  hhoCell         : hho Cell
! In hhoData          : information about the HHO formulation
! Out gradfull        : full gradient
! Out stab            : stabilization
! --------------------------------------------------------------------------------------------------
!
! --- Local variables
!
        integer :: cbs, fbs, total_dofs, gbs
!
        call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)
!
! -------- Reload gradient
        call gradfull%initialize(gbs, total_dofs, 0.d0)
        call gradfull%read('PCHHOGT', ASTER_FALSE)
!
! -------- Reload stabilization
        if (present(stab)) then
            call stab%initialize(total_dofs, total_dofs, 0.d0)
            call stab%read('PCHHOST', ASTER_TRUE)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalcOpTher(hhoCell, hhoData, gradfull, stab)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(out) :: gradfull
        type(HHO_matrix), intent(out), optional :: stab
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Thermics
!
! Compute operators for thermics
!
! --------------------------------------------------------------------------------------------------
!
! In  hhoCell         : hho Cell
! In hhoData          : information about the HHO formulation
! Out gradfull        : full gradient
! Out stab            : stabilization
! --------------------------------------------------------------------------------------------------
!
        type(HHO_matrix) :: gradrec_scal
!
! ----- Compute Gradient reconstruction
        call hhoGradRecFullVec(hhoCell, hhoData, gradfull)
!
! ----- Compute Stabilizatiion
        if (present(stab)) then
            if (hhoData%cell_degree() <= hhoData%face_degree()) then
                call hhoGradRecVec(hhoCell, hhoData, gradrec_scal)
                call hhoStabScal(hhoCell, hhoData, gradrec_scal, stab)
                call gradrec_scal%free()
            else if (hhoData%cell_degree() == (hhoData%face_degree()+1)) then
                call hdgStabScal(hhoCell, hhoData, stab)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeRhsRigiTher(hhoCell, stress, weight, BSCEval, gbs, bT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: stress(3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer, intent(in) :: gbs
        real(kind=8), intent(inout) :: bT(MSIZE_CELL_VEC)
!
! ------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the scalar product bT += (stress, gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In stress       : stress tensor
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs          : number of rows of bT
!   Out bT          : contribution of bt
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_stress(3)
        integer :: i, gbs_cmp, deca
! --------------------------------------------------------------------------------------------------
!
        gbs_cmp = gbs/hhoCell%ndim
        qp_stress = weight*stress
! -------- Compute scalar_product of (stress, gphi)_T
! -------- (RAPPEL: the composents of the gradient are stored by rows)
        deca = 0
        do i = 1, hhoCell%ndim
            call daxpy_1(gbs_cmp, qp_stress(i), BSCEval, bT(deca+1))
            deca = deca+gbs_cmp
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeRhsMassTher(beta, weight, BSCEval, cbs, rhs)
!
        implicit none
!
        real(kind=8), intent(in) :: beta
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer, intent(in) :: cbs
        real(kind=8), intent(inout) :: rhs(MSIZE_CELL_SCAL)
!
! ------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the scalar product rhs += (beta, vT)_T at a quadrature point
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In cbs          : number of rows of bT
!   Out bT          : contribution of bt
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_temp
! --------------------------------------------------------------------------------------------------
!
        qp_temp = weight*beta
        call daxpy_1(cbs, qp_temp, BSCEval, rhs)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsRigiTher(hhoCell, module_tang, weight, BSCEval, gbs, &
                                     AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: module_tang(3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer, intent(in) :: gbs
        real(kind=8), intent(inout) :: AT(MSIZE_CELL_VEC, MSIZE_CELL_VEC)
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the scalar product AT += (module_tang:gphi, gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs          : number of rows of bT
!   Out AT          : contribution of bt
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Agphi(MSIZE_CELL_VEC, 3)
        integer :: i, k, gbs_cmp, col
! --------------------------------------------------------------------------------------------------
!
        gbs_cmp = gbs/hhoCell%ndim
! --------- Eval (A : gphi)_T
        call hhoComputeAgphi(hhoCell, module_tang, BSCEval, gbs, weight, &
                             qp_Agphi)
!
! -------- Compute scalar_product of (A_gphi, gphi)_T
        col = 1
        do i = 1, hhoCell%ndim
            do k = 1, gbs_cmp
                call daxpy_1(gbs, BSCEval(k), qp_Agphi(:, i), AT(:, col))
                col = col+1
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsMassTher(cp, weight, BSCEval, cbs, lhs)
!
        implicit none
!
        real(kind=8), intent(in) :: cp
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        integer, intent(in) :: cbs
        type(HHO_matrix), intent(inout) :: lhs
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the scalar product lhs += (cp * cphi, cphi)_T at a quadrature point
!   In cp           : rho * Cp
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In cbs          : number of rows of bT
!   Out rhs         : contribution of bt
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: coeff
        blas_int :: b_incx, b_lda, b_n
! --------------------------------------------------------------------------------------------------
!
        coeff = cp*weight
        b_n = to_blas_int(cbs)
        b_incx = to_blas_int(1)
        b_lda = to_blas_int(lhs%max_nrows)
        call dsyr('U', b_n, coeff, BSCEval, b_incx, lhs%m, b_lda)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeAgphi(hhoCell, module_tang, BSCEval, gbs, weight, Agphi)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        integer, intent(in) :: gbs
        real(kind=8), intent(in) :: module_tang(3, 3)
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(out) :: Agphi(MSIZE_CELL_VEC, 3)
!
! -----------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product (module_tang, gphi)_T
!
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto_plastic moduli
!
!   In BSCEval      : Basis of one composant gphi
!   In gbs          : number of cols of Aphi
!   In weight       : quadrature weight
!   Out Agphi       : matriw of scalar product
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_module_tang(3, 3), qp_mod_vec(3)
        integer :: i, row, gbs_cmp, dim
        blas_int :: b_incx, b_incy, b_lda, b_m, b_n
! --------------------------------------------------------------------------------------------------
!
        Agphi = 0.d0
        dim = hhoCell%ndim
        gbs_cmp = gbs/dim
        qp_module_tang = weight*module_tang
!
        row = 1
        do i = 1, dim
            qp_mod_vec = qp_module_tang(:, i)
            b_lda = to_blas_int(gbs_cmp)
            b_m = to_blas_int(gbs_cmp)
            b_n = to_blas_int(dim)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dger(b_m, b_n, 1.d0, BSCEval, b_incx, &
                      qp_mod_vec, b_incy, Agphi(row:(row+gbs_cmp-1), 1:dim), b_lda)
            row = row+gbs_cmp
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function LambdaMax(fami, imate, npg, ndim, time, &
                       temp_eval) result(coeff)
!
        implicit none
!
        character(len=8), intent(in) :: fami
        integer, intent(in) :: imate, npg, ndim
        real(kind=8), intent(in) :: time, temp_eval
        real(kind=8) :: coeff
!
! --------------------------------------------------------------------------------------------------
!
!   Compute the average Young modulus
!   In fami         : familly of quadrature points (of hhoQuadCellRigi)
!   In npg          : number of quadrature points
!   In imate        : materiau code
! --------------------------------------------------------------------------------------------------
!
        character(len=16) :: nomres(3)
        real(kind=8) :: valres(3)
        integer :: icodre(3), ipg
        integer, parameter :: spt = 1
        character(len=32) :: phenom
        real(kind=8) :: lambda, lambor(3)
!
        coeff = 1.d0
        call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
!
        do ipg = 1, npg
            if (phenom .eq. 'THER') then
                nomres(1) = 'LAMBDA'
                call rcvalb(fami, ipg, spt, '+', zi(imate), &
                            ' ', phenom, 1, 'INST', [time], &
                            1, nomres, valres, icodre, 1)
                lambda = valres(1)
                coeff = coeff+lambda
            else if (phenom .eq. 'THER_NL') then
                nomres(1) = 'LAMBDA'
                call rcvalb(fami, ipg, spt, '+', zi(imate), &
                            ' ', phenom, 1, 'TEMP', [temp_eval], &
                            1, nomres, valres, icodre, 1)
                lambda = valres(1)
                coeff = coeff+lambda
            else if (phenom .eq. 'THER_ORTH') then
                nomres(1) = 'LAMBDA_L'
                nomres(2) = 'LAMBDA_T'
                nomres(3) = 'LAMBDA_N'
                call rcvalb(fami, ipg, spt, '+', zi(imate), &
                            ' ', phenom, 1, 'INST', [time], &
                            ndim, nomres, valres, icodre, 1)
                lambor(1) = valres(1)
                lambor(2) = valres(2)
                lambor(3) = 0.d0
                if (ndim == 3) lambor(3) = valres(3)
                coeff = coeff+maxval(lambor)
            else if (phenom .eq. 'THER_NL_ORTH') then
                nomres(1) = 'LAMBDA_L'
                nomres(2) = 'LAMBDA_T'
                nomres(3) = 'LAMBDA_N'
                call rcvalb(fami, ipg, spt, '+', zi(imate), &
                            ' ', phenom, 1, 'TEMP', [temp_eval], &
                            ndim, nomres, valres, icodre, 1)
                lambor(1) = valres(1)
                lambor(2) = valres(2)
                lambor(3) = 0.d0
                if (ndim == 3) lambor(3) = valres(3)
                coeff = coeff+maxval(lambor)
            else
                call utmess('F', 'ELEMENTS2_63')
            end if
!
        end do
!
        coeff = coeff/real(npg, kind=8)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalcStabCoeffTher(hhoData, fami, nbQuadPoints, ndim, time)
!
        implicit none
!
        type(HHO_Data), intent(inout) :: hhoData
        character(len=8) :: fami
        integer, intent(in) :: nbQuadPoints, ndim
        real(kind=8), intent(in) :: time
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - Evaluate stabilzation coefficient
!
! In hhoData          : information about the HHO formulation
! --------------------------------------------------------------------------------------------------
!
! --- Local variables
!
        integer :: imate
!
        if (hhoData%adapt()) then
            call jevech('PMATERC', 'L', imate)
            call hhoData%setCoeffStab(10.d0*LambdaMax(fami, imate, nbQuadPoints, ndim, time, 20.0 &
                                                      ))
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeBehaviourTher(phenom, fami, kpg, ndim, time, &
                                       imate, coorpg, temp, dtemp, fluxglo, &
                                       Kglo)
!
        implicit none
!
        character(len=32), intent(in) :: phenom
        character(len=8), intent(in) :: fami
        integer, intent(in) :: imate, ndim, kpg
        real(kind=8), intent(in) :: time, dtemp(3), temp, coorpg(3)
        real(kind=8), intent(out) :: fluxglo(3), Kglo(3, 3)
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - Integrate behaviour - Fluxes and derivative
! --------------------------------------------------------------------------------------------------
!
        call nlcomp(phenom, fami, kpg, imate, ndim, &
                    coorpg, time, temp, Kglo, dtp_=dtemp, &
                    fluglo_=fluxglo)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeRhoCpTher(phenom, fami, poum, kpg, time, &
                                   imate, cp)
!
        implicit none
!
        character(len=32), intent(in) :: phenom
        character(len=8), intent(in) :: fami, poum
        integer, intent(in) :: imate, kpg
        real(kind=8), intent(in) :: time
        real(kind=8), intent(out) :: cp
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - Compute rho_cp
!
! In hhoData          : information about the HHO formulation
! --------------------------------------------------------------------------------------------------
!
! --- Local variables
!
        integer :: icod(1)
        integer, parameter :: spt = 1
        real(kind=8) :: cp_(1)
!
! --- Eval rho_cp
!
        if (phenom .eq. 'THER') then
            call rcvalb(fami, kpg, spt, poum, zi(imate), &
                        ' ', phenom, 1, 'INST', [time], &
                        1, 'RHO_CP', cp_, icod(1), 1)
        else if (phenom .eq. 'THER_ORTH') then
            call rcvalb(fami, kpg, spt, poum, zi(imate), &
                        ' ', phenom, 1, 'INST', [time], &
                        1, 'RHO_CP', cp_, icod(1), 1)
        else
            call utmess('F', 'ELEMENTS2_63')
        end if
!
        cp = cp_(1)
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeBetaTher(comp, ifon, temp, beta, dbeta)
!
        implicit none
!
        character(len=16), intent(in) :: comp
        integer, intent(in) :: ifon
        real(kind=8), intent(in) :: temp
        real(kind=8), intent(out) :: beta, dbeta
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - Compute beta and d_beta
!
! --------------------------------------------------------------------------------------------------
!
! --- Eval beta and derivative
!
        if (comp(1:5) .eq. 'THER_') then
            call rcfode(ifon, temp, beta, dbeta)
            if (comp(1:9) .eq. 'THER_HYDR') then
                ASSERT(ASTER_FALSE)
            end if
        else if (comp(1:5) .eq. 'SECH_') then
            beta = temp
            dbeta = 1.d0
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
end module
