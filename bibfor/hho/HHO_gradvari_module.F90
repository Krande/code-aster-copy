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
!
! aslint: disable=C9992

! WARNING: Some big arrays are larger than limit set by '-fmax-stack-var-size='.
! The 'save' attribute has been added. They *MUST NOT* been accessed concurrently.

module HHO_GV_module
!
    use NonLin_Datastructure_type
    use Behaviour_type
    use Behaviour_module
    use HHO_basis_module
    use HHO_compor_module
    use HHO_Dirichlet_module
    use HHO_eval_module
    use HHO_LargeStrainMeca_module
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_SmallStrainMeca_module
    use HHO_stabilization_module, only: hhoStabVec, hdgStabVec, hhoStabSymVec
    use HHO_Ther_module
    use HHO_type
    use HHO_utils_module
    use HHO_gradrec_module, only: hhoGradRecVec, hhoGradRecFullMat, hhoGradRecSymFullMat
    use HHO_gradrec_module, only: hhoGradRecSymMat, hhoGradRecFullMatFromVec
!
    implicit none
!
    private
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/desymt46.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/lagmodtonommod.h"
#include "asterfort/lcdetf.h"
#include "asterfort/nmcomp.h"
#include "asterfort/sigtopk1.h"
#include "asterfort/poslog.h"
#include "asterfort/prelog.h"
#include "asterfort/readVector.h"
#include "asterfort/rcvalb.h"
#include "asterfort/deflg4.h"
#include "asterfort/prodmt.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dsymv.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
#include "blas/dsyr.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - mechanics - GRAD_VARI
!
! Specific routines for mechanics
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: HHO_GV_State
    public :: hhoGradVariLC, hhoCalcOpGv
    private :: check_behavior, gdef_log, hhoAssGVRhs, hhoAssGVLhs, numGVMap
    private :: initialize_gv, hhoCalcStabCoeffGV
!
! --------------------------------------------------------------------------------------------------
!
    type HHO_GV_State
!
        aster_logical :: l_debug = ASTER_FALSE
! ----- GRAD_VARI
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: vari_prev = 0.d0
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: vari_curr = 0.d0
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: vari_incr = 0.d0
! ----- LAGR
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: lagv_prev = 0.d0
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: lagv_curr = 0.d0
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: lagv_incr = 0.d0
! ----- Gradient reconstruction and stabilisation
        real(kind=8) :: grad(MSIZE_CELL_VEC, MSIZE_TDOFS_SCAL) = 0.d0
        real(kind=8) :: stab(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL) = 0.d0
! ----- member function
    contains
        procedure, pass :: initialize => initialize_gv
    end type HHO_GV_State
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGradVariLC(hhoCell, hhoData, hhoQuadCellRigi, hhoMecaState, hhoComporState, &
                             hhoGVState, lhs, rhs)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuadCellRigi
        type(HHO_Meca_State), intent(in) :: hhoMecaState
        type(HHO_Compor_State), intent(inout) :: hhoComporState
        type(HHO_GV_State), intent(in) :: hhoGVState
        real(kind=8), intent(out) :: lhs(MSIZE_TDOFS_MIX, MSIZE_TDOFS_MIX)
        real(kind=8), intent(out) :: rhs(MSIZE_TDOFS_MIX)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the local contribution for GRAD_VARI
!   In hhoCell      : the current HHO Cell
!   In hhoData       : information on HHO methods
!   In hhoQuadCellRigi : quadrature rules from the rigidity family
!   Out lhs         : local contribution (lhs)
!   Out rhs         : local contribution (rhs)
!
!   phi = phi^e(eps) + phi^p(p) + c/2(grad vari, grad_vari) + phi^nl(vari, lagv, p)
!   eps_gene = (eps, vari, lagv, grad vari)
!   sig_gene = (sig, dphi^nl/dvari, dphi^nl/dlagv, c grad vari)
!
!   Attention, il faut que la quadrature soit suffisante pour tout les termes comme
!   (c_phi, c_phi), (c_phi, g_phi), (g_phi, g_phi)
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        type(Behaviour_Integ) :: BEHinteg
        real(kind=8), dimension(MSIZE_CELL_MAT) :: mk_bT, G_prev_coeff, G_curr_coeff
        real(kind=8), dimension(MSIZE_CELL_VEC) :: gv_bT, GV_prev_coeff, GV_curr_coeff
        real(kind=8) :: G_prev(3, 3), G_curr(3, 3)
        real(kind=8) :: GV_prev(3), GV_curr(3)
        real(kind=8) :: F_prev(3, 3), F_curr(3, 3), Pk1(3, 3)
        real(kind=8) :: var_prev, var_curr, lag_prev, lag_curr
        real(kind=8) :: sig_vari, sig_lagv, sig_gv(3), dsv_dv, dsv_dl, dsl_dl, dsgv_dgv(3, 3)
        real(kind=8) :: dsv_dF(3, 3), dsl_dF(3, 3), dsv_dEps(6), dsl_dEps(6)
        real(kind=8) :: Eps_prev(6), Eps_curr(6), Cauchy(6)
        real(kind=8) :: dPK1_dF(3, 3, 3, 3), dPK1_dv(3, 3), dPK1_dl(3, 3)
        real(kind=8) :: dSig_dEps(6, 6), dSig_dv(6), dSig_dl(6)
        real(kind=8) :: jac_prev, jac_curr, coorpg(3), weight, coeff, mk_stab, gv_stab
        real(kind=8) :: BSCEvalG(MSIZE_CELL_SCAL), BSCEval(MSIZE_CELL_SCAL)
        real(kind=8) :: mk_AT(MSIZE_CELL_MAT, MSIZE_CELL_MAT)
        real(kind=8) :: mk_TMP(MSIZE_CELL_MAT, MSIZE_TDOFS_VEC)
        real(kind=8) :: gv_AT(MSIZE_CELL_VEC, MSIZE_CELL_VEC)
        real(kind=8) :: gv_TMP(MSIZE_CELL_VEC, MSIZE_TDOFS_SCAL)
        real(kind=8) :: mv_AT(MSIZE_CELL_MAT, MSIZE_CELL_SCAL)
        real(kind=8) :: ml_AT(MSIZE_CELL_MAT, MSIZE_CELL_SCAL)
        real(kind=8) :: vm_AT(MSIZE_CELL_SCAL, MSIZE_CELL_MAT)
        real(kind=8) :: lm_AT(MSIZE_CELL_SCAL, MSIZE_CELL_MAT)
        real(kind=8) :: rhs_vari(MSIZE_TDOFS_SCAL), rhs_lagv(MSIZE_CELL_SCAL)
        real(kind=8) :: rhs_mk(MSIZE_TDOFS_VEC)
        real(kind=8) :: lhs_mm(MSIZE_TDOFS_VEC, MSIZE_TDOFS_VEC)
        real(kind=8) :: lhs_mv(MSIZE_TDOFS_VEC, MSIZE_TDOFS_SCAL)
        real(kind=8) :: lhs_ml(MSIZE_TDOFS_VEC, MSIZE_CELL_SCAL)
        real(kind=8) :: lhs_vm(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_VEC)
        real(kind=8) :: lhs_vv(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL)
        real(kind=8) :: lhs_vl(MSIZE_TDOFS_SCAL, MSIZE_CELL_SCAL)
        real(kind=8) :: lhs_lm(MSIZE_CELL_SCAL, MSIZE_TDOFS_VEC)
        real(kind=8) :: lhs_lv(MSIZE_CELL_SCAL, MSIZE_TDOFS_SCAL)
        real(kind=8) :: lhs_ll(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
        integer :: mapMeca(MSIZE_TDOFS_VEC), mapVari(MSIZE_TDOFS_SCAL), mapLagv(MSIZE_CELL_SCAL)
        integer :: mk_cbs, mk_fbs, mk_total_dofs, mk_gbs, mk_gbs_sym, mk_gbs_cmp
        integer :: gv_cbs, gv_fbs, gv_total_dofs, gv_gbs
        integer :: cod(27), ipg, j, mk_gbs_tot
        aster_logical :: l_lhs, l_rhs, forc_noda
        blas_int :: b_incx, b_incy, b_n
        blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m
! --------------------------------------------------------------------------------------------------
!
        cod = 0
        lhs = 0.d0
        rhs = 0.d0
!
! ------ number of dofs
!
        call hhoMecaNLDofs(hhoCell, hhoData, mk_cbs, mk_fbs, mk_total_dofs, &
                           mk_gbs, mk_gbs_sym)
        call hhoTherNLDofs(hhoCell, hhoData, gv_cbs, gv_fbs, gv_total_dofs, &
                           gv_gbs)
        if (hhoComporState%l_largestrain) then
            mk_gbs_tot = mk_gbs
        else
            mk_gbs_tot = mk_gbs_sym
        end if
        mk_gbs_cmp = gv_gbs
!
! ------ initialization
!
        mk_bT = 0.d0
        mk_AT = 0.d0
        mk_TMP = 0.d0
        G_prev_coeff = 0.d0
        G_curr_coeff = 0.d0
        rhs_mk = 0.d0
        lhs_mm = 0.d0
        lhs_mv = 0.d0
        lhs_ml = 0.d0
!
        gv_bT = 0.d0
        gv_AT = 0.d0
        gv_tmp = 0.d0
        GV_prev_coeff = 0.d0
        GV_curr_coeff = 0.d0
        rhs_vari = 0.d0
        rhs_lagv = 0.d0
        lhs_vm = 0.d0
        lhs_vv = 0.d0
        lhs_vl = 0.d0
        lhs_lm = 0.d0
        lhs_lv = 0.d0
        lhs_ll = 0.d0
        mv_AT = 0.d0
        ml_AT = 0.d0
        vm_AT = 0.d0
        lm_AT = 0.d0
!
! ----- Type of behavior
!
        call check_behavior(hhoComporState)
!
! ----- Initialisation of behaviour datastructure
!
        call behaviourInit(BEHinteg)
!
        forc_noda = hhoComporState%option == "FORC_NODA"
        l_lhs = L_MATR(hhoComporState%option)
        l_rhs = L_VECT(hhoComporState%option) .or. forc_noda
!
! ----- init basis
!
        call hhoBasisCell%initialize(hhoCell)
!
! ----- compute G_prev = gradrec * depl_prev (sym or not)
!
        b_lda = to_blas_int(MSIZE_CELL_MAT)
        b_m = to_blas_int(mk_gbs_tot)
        b_n = to_blas_int(mk_total_dofs)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, 1.d0, hhoMecaState%grad, &
                   b_lda, hhoMecaState%depl_prev, b_incx, 0.d0, G_prev_coeff, &
                   b_incy)
!
! ----- compute G_curr = gradrec * depl_curr (sym or not)
!
        b_lda = to_blas_int(MSIZE_CELL_MAT)
        b_m = to_blas_int(mk_gbs_tot)
        b_n = to_blas_int(mk_total_dofs)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, 1.d0, hhoMecaState%grad, &
                   b_lda, hhoMecaState%depl_curr, b_incx, 0.d0, G_curr_coeff, &
                   b_incy)
!
! ----- compute GV_prev = gradrec * vari_prev
!
        b_lda = to_blas_int(MSIZE_CELL_VEC)
        b_m = to_blas_int(gv_gbs)
        b_n = to_blas_int(gv_total_dofs)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, 1.d0, hhoGVState%grad, &
                   b_lda, hhoGVState%vari_prev, b_incx, 0.d0, GV_prev_coeff, &
                   b_incy)
!
! ----- compute GV_curr = gradrec * vari_curr
!
        b_lda = to_blas_int(MSIZE_CELL_VEC)
        b_m = to_blas_int(gv_gbs)
        b_n = to_blas_int(gv_total_dofs)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dgemv('N', b_m, b_n, 1.d0, hhoGVState%grad, &
                   b_lda, hhoGVState%vari_curr, b_incx, 0.d0, GV_curr_coeff, &
                   b_incy)
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCellRigi%nbQuadPoints
            coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
            weight = hhoQuadCellRigi%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
            call hhoBasisCell%BSEval(hhoCell, coorpg(1:3), 0, hhoData%grad_degree(), BSCEvalG)
            call hhoBasisCell%BSEval(hhoCell, coorpg(1:3), 0, hhoData%cell_degree(), BSCEval)
!
! --------- Eval gradient at T- and T+
!
            if (hhoComporState%l_largestrain) then
                G_prev = hhoEvalMatCell( &
                         hhoCell, hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), G_prev_coeff, &
                         mk_gbs &
                         )
                G_curr = hhoEvalMatCell( &
                         hhoCell, hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), G_curr_coeff, &
                         mk_gbs &
                         )
!
! --------- Eval gradient of the deformation at T- and T+
!
                call hhoCalculF(hhoCell%ndim, G_prev, F_prev)
                call hhoCalculF(hhoCell%ndim, G_curr, F_curr)
!
! -------- Check the jacobian jac >= r8prem
! -------- be carrefull with c_plan, I don't know the result
!
                call lcdetf(hhoCell%ndim, F_prev, jac_prev)
                cod(ipg) = merge(1, 0, jac_prev .le. r8prem())
                if (cod(ipg) .ne. 0) goto 999
!
                call lcdetf(hhoCell%ndim, F_curr, jac_curr)
                cod(ipg) = merge(1, 0, jac_curr .le. r8prem())
                if (cod(ipg) .ne. 0) goto 999
            else
                Eps_prev = hhoEvalSymMatCell( &
                           hhoCell, hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), &
                           G_prev_coeff, mk_gbs_sym &
                           )
                Eps_curr = hhoEvalSymMatCell( &
                           hhoCell, hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), &
                           G_curr_coeff, mk_gbs_sym &
                           )
            end if
!
            GV_prev = hhoEvalVecCell( &
                      hhoCell, hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), GV_prev_coeff, &
                      gv_gbs &
                      )
            GV_curr = hhoEvalVecCell( &
                      hhoCell, hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), GV_curr_coeff, &
                      gv_gbs &
                      )
!
            var_prev = hhoEvalScalCell( &
                       hhoCell, hhoBasisCell, hhoData%cell_degree(), coorpg(1:3), &
                       hhoGVState%vari_prev, gv_cbs &
                       )
            var_curr = hhoEvalScalCell( &
                       hhoCell, hhoBasisCell, hhoData%cell_degree(), coorpg(1:3), &
                       hhoGVState%vari_curr, gv_cbs &
                       )
!
            lag_prev = hhoEvalScalCell( &
                       hhoCell, hhoBasisCell, hhoData%cell_degree(), coorpg(1:3), &
                       hhoGVState%lagv_prev, gv_cbs &
                       )
            lag_curr = hhoEvalScalCell( &
                       hhoCell, hhoBasisCell, hhoData%cell_degree(), coorpg(1:3), &
                       hhoGVState%lagv_curr, gv_cbs &
                       )
!
! ------- Compute behavior
!
            if (forc_noda) then
                call forc_noda_stress(hhoComporState, hhoCell%ndim, ipg, F_curr, Pk1, &
                                      Cauchy, sig_vari, sig_lagv, sig_gv)
            else if (hhoComporState%l_largestrain) then
                call gdef_log(BEHinteg, hhoComporState, hhoCell%ndim, ipg, &
                              hhoMecaState%time_prev, hhoMecaState%time_curr, F_prev, F_curr, &
                              var_prev, var_curr, lag_prev, lag_curr, GV_prev, &
                              GV_curr, Pk1, sig_vari, sig_lagv, sig_gv, &
                              dPK1_dF, dPK1_dv, dPK1_dl, dsv_dF, dsv_dv, &
                              dsv_dl, dsl_dF, dsl_dl, dsgv_dgv, cod(ipg))
            else
                call petit(BEHinteg, hhoComporState, hhoCell%ndim, ipg, hhoMecaState%time_prev, &
                           hhoMecaState%time_curr, Eps_prev, Eps_curr, var_prev, var_curr, &
                           lag_prev, lag_curr, GV_prev, GV_curr, Cauchy, &
                           sig_vari, sig_lagv, sig_gv, dSig_dEps, dSig_dv, &
                           dSig_dl, dsv_dEps, dsv_dv, dsv_dl, dsl_dEps, &
                           dsl_dl, dsgv_dgv, cod(ipg))
            end if
!
! -------- Test the code of the LDC
!
            if (cod(ipg) .eq. 1) goto 999
!
! ------- Compute rhs
!
            if (l_rhs) then
                if (hhoComporState%l_largestrain) then
! ---------- += weight * (PK1, g_phi)
                    call hhoComputeRhsLarge(hhoCell, Pk1, weight, BSCEvalG, mk_gbs, &
                                            mk_bT)
                else
! ---------- += weight * (Cauchy, gs_phi)
                    call hhoComputeRhsSmall(hhoCell, Cauchy, weight, BSCEvalG, mk_gbs_cmp, &
                                            mk_bT)
                end if
! ---------- += weight * (sig_gv, g_phi)
                call hhoComputeRhsRigiTher(hhoCell, sig_gv, weight, BSCEvalG, gv_gbs, &
                                           gv_bT)
! ---------- += weight * (sig_vari, c_phi)
                call hhoComputeRhsMassTher(1.d0, sig_vari, weight, BSCEval, gv_cbs, &
                                           rhs_vari)
! ---------- += weight * (sig_lagv, c_phi)
                call hhoComputeRhsMassTher(1.d0, sig_lagv, weight, BSCEval, gv_cbs, &
                                           rhs_lagv)
            end if
!
! ------- Compute lhs
!
            if (l_lhs) then
                if (hhoComporState%l_largestrain) then
! ---------- += weight * (dPK1_dF : g_phi, g_phi)
                    call hhoComputeLhsLarge(hhoCell, dPK1_dF, weight, BSCEvalG, mk_gbs, &
                                            mk_AT)
! ---------- += weight * (g_phi, dPK1_dv : c_phi) -> lhs_mv
                    call hhoComputeLhsLargeMG(hhoCell, dPK1_dv, weight, BSCEval, gv_cbs, &
                                              BSCEvalG, mk_gbs, mv_AT)
! ---------- += weight * (g_phi, dPK1_dl : c_phi) -> lhs_ml
                    call hhoComputeLhsLargeMG(hhoCell, dPK1_dl, weight, BSCEval, gv_cbs, &
                                              BSCEvalG, mk_gbs, ml_AT)
! ---------- += weight * (dsv_dF : g_phi, c_phi) -> lhs_vm
                    call hhoComputeLhsLargeGM(hhoCell, dsv_dF, weight, BSCEval, gv_cbs, &
                                              BSCEvalG, mk_gbs, vm_AT)
! ---------- += weight * (dsl_dF : g_phi, c_phi) -> lhs_lm
                    call hhoComputeLhsLargeGM(hhoCell, dsl_dF, weight, BSCEval, gv_cbs, &
                                              BSCEvalG, mk_gbs, lm_AT)
                else
! ---------- += weight * (dSig_deps : gs_phi, gs_phi)
                    call hhoComputeLhsSmall(hhoCell, dSig_deps, weight, BSCEvalG, mk_gbs_sym, &
                                            mk_gbs_cmp, mk_AT)
! ---------- TODO: Add missing terms
! ---------- += weight * (gs_phi, dSig_dv : c_phi) -> lhs_mv
!                     call hhoComputeLhsSmallMG(hhoCell, dSig_dv, weight, BSCEval, gv_cbs, &
!                                             BSCEvalG, mk_gbs, mv_AT)
! ! ---------- += weight * (gs_phi, dSig_dl : c_phi) -> lhs_ml
!                     call hhoComputeLhsSmallMG(hhoCell, dSig_dl, weight, BSCEval, gv_cbs, &
!                                             BSCEvalG, mk_gbs, ml_AT)
! ! ---------- += weight * (dsv_dEps : gs_phi, c_phi) -> lhs_vm
!                     call hhoComputeLhsSmallGM(hhoCell, dsv_dEps, weight, BSCEval, gv_cbs, &
!                                             BSCEvalG, mk_gbs, vm_AT)
! ! ---------- += weight * (dsl_dEps : g_phi, c_phi) -> lhs_lm
!                     call hhoComputeLhsSmallGM(hhoCell, dsl_dEps, weight, BSCEval, gv_cbs, &
!                                             BSCEvalG, mk_gbs, lm_AT)
                end if
! ---------- += weight * (dgv_dv : g_phi, g_phi)
                call hhoComputeLhsRigiTher(hhoCell, dsgv_dgv, weight, BSCEvalG, gv_gbs, &
                                           gv_AT)
! ---------- += weight * (dsv_dv : c_phi, c_phi)
                coeff = weight*dsv_dv
                b_n = to_blas_int(gv_cbs)
                b_incx = to_blas_int(1)
                b_lda = to_blas_int(MSIZE_TDOFS_SCAL)
                call dsyr('U', b_n, coeff, BSCEval, b_incx, &
                          lhs_vv, b_lda)
! ---------- += weight * (dsv_dl : c_phi, c_phi)
                coeff = weight*dsv_dl
                b_n = to_blas_int(gv_cbs)
                b_incx = to_blas_int(1)
                b_lda = to_blas_int(MSIZE_TDOFS_SCAL)
                call dsyr('U', b_n, coeff, BSCEval, b_incx, &
                          lhs_vl, b_lda)
! ---------- += weight * (dsl_dl : c_phi, c_phi)
                call hhoComputeLhsMassTher(dsl_dl, weight, BSCEval, gv_cbs, lhs_ll)
            end if
!
        end do
!
        call numGVMap(hhoCell, hhoData, mapMeca, mapVari, mapLagv)
        call hhoCalcStabCoeffMeca(hhoData, hhoComporState%fami, hhoMecaState%time_curr, &
                                  hhoQuadCellRigi)
        mk_stab = hhoData%coeff_stab()
        gv_stab = hhoCalcStabCoeffGV(hhoComporState%fami, hhoQuadCellRigi%nbQuadPoints)
!
! ------- Compute rhs
!
        if (l_rhs) then
! ----- compute rhs += Gradrec**T * bT
            b_lda = to_blas_int(MSIZE_CELL_MAT)
            b_m = to_blas_int(mk_gbs_tot)
            b_n = to_blas_int(mk_total_dofs)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dgemv('T', b_m, b_n, 1.d0, hhoMecaState%grad, &
                       b_lda, mk_bT, b_incx, 1.d0, rhs_mk, &
                       b_incy)
! ----- compute rhs += stab
            b_lda = to_blas_int(MSIZE_TDOFS_VEC)
            b_n = to_blas_int(mk_total_dofs)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dsymv('U', b_n, mk_stab, hhoMecaState%stab, b_lda, &
                       hhoMecaState%depl_curr, b_incx, 1.d0, rhs_mk, b_incy)
! ----- compute rhs += Gradrec**T * bT
            b_lda = to_blas_int(MSIZE_CELL_VEC)
            b_m = to_blas_int(gv_gbs)
            b_n = to_blas_int(gv_total_dofs)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dgemv('T', b_m, b_n, 1.d0, hhoGVState%grad, &
                       b_lda, gv_bT, b_incx, 1.d0, rhs_vari, &
                       b_incy)
! ----- compute rhs += stab
            b_lda = to_blas_int(MSIZE_TDOFS_SCAL)
            b_n = to_blas_int(gv_total_dofs)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dsymv('U', b_n, gv_stab, hhoGVState%stab, b_lda, &
                       hhoGVState%vari_curr, b_incx, 1.d0, rhs_vari, b_incy)
! ----- assembly
            call hhoAssGVRhs(hhoCell, hhoData, mapMeca, mapVari, mapLagv, &
                             rhs_mk, rhs_vari, rhs_lagv, rhs)
        end if
!
! ------- Compute lhs
!
        if (l_lhs) then
! ----- Add symmetry
!
            call hhoCopySymPartMat('U', lhs_vv(1:gv_cbs, 1:gv_cbs))
            call hhoCopySymPartMat('U', lhs_vl(1:gv_cbs, 1:gv_cbs))
            call hhoCopySymPartMat('U', lhs_ll(1:gv_cbs, 1:gv_cbs))
!
! ----- Add gradient: += gradrec**T * AT * gradrec
! ----- step1: TMP = AT * gradrec
            b_ldc = to_blas_int(MSIZE_CELL_MAT)
            b_ldb = to_blas_int(MSIZE_CELL_MAT)
            b_lda = to_blas_int(MSIZE_CELL_MAT)
            b_m = to_blas_int(mk_gbs_tot)
            b_n = to_blas_int(mk_total_dofs)
            b_k = to_blas_int(mk_gbs_tot)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       1.d0, mk_AT, b_lda, hhoMecaState%grad, b_ldb, &
                       0.d0, mk_TMP, b_ldc)
            b_ldc = to_blas_int(MSIZE_CELL_VEC)
            b_ldb = to_blas_int(MSIZE_CELL_VEC)
            b_lda = to_blas_int(MSIZE_CELL_VEC)
            b_m = to_blas_int(gv_gbs)
            b_n = to_blas_int(gv_total_dofs)
            b_k = to_blas_int(gv_gbs)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       1.d0, gv_AT, b_lda, hhoGVState%grad, b_ldb, &
                       0.d0, gv_TMP, b_ldc)
! ----- step2: lhs += gradrec**T * TMP
            b_ldc = to_blas_int(MSIZE_TDOFS_VEC)
            b_ldb = to_blas_int(MSIZE_CELL_MAT)
            b_lda = to_blas_int(MSIZE_CELL_MAT)
            b_m = to_blas_int(mk_total_dofs)
            b_n = to_blas_int(mk_total_dofs)
            b_k = to_blas_int(mk_gbs_tot)
            call dgemm('T', 'N', b_m, b_n, b_k, &
                       1.d0, hhoMecaState%grad, b_lda, mk_TMP, b_ldb, &
                       0.d0, lhs_mm, b_ldc)
            b_ldc = to_blas_int(MSIZE_TDOFS_SCAL)
            b_ldb = to_blas_int(MSIZE_CELL_VEC)
            b_lda = to_blas_int(MSIZE_CELL_VEC)
            b_m = to_blas_int(gv_total_dofs)
            b_n = to_blas_int(gv_total_dofs)
            b_k = to_blas_int(gv_gbs)
            call dgemm('T', 'N', b_m, b_n, b_k, &
                       1.d0, hhoGVState%grad, b_lda, gv_TMP, b_ldb, &
                       1.d0, lhs_vv, b_ldc)
!
            b_ldc = to_blas_int(MSIZE_TDOFS_VEC)
            b_ldb = to_blas_int(MSIZE_CELL_MAT)
            b_lda = to_blas_int(MSIZE_CELL_MAT)
            b_m = to_blas_int(mk_total_dofs)
            b_n = to_blas_int(gv_cbs)
            b_k = to_blas_int(mk_gbs_tot)
            call dgemm('T', 'N', b_m, b_n, b_k, &
                       1.d0, hhoMecaState%grad, b_lda, mv_AT, b_ldb, &
                       0.d0, lhs_mv, b_ldc)
            b_ldc = to_blas_int(MSIZE_TDOFS_VEC)
            b_ldb = to_blas_int(MSIZE_CELL_MAT)
            b_lda = to_blas_int(MSIZE_CELL_MAT)
            b_m = to_blas_int(mk_total_dofs)
            b_n = to_blas_int(gv_cbs)
            b_k = to_blas_int(mk_gbs_tot)
            call dgemm('T', 'N', b_m, b_n, b_k, &
                       1.d0, hhoMecaState%grad, b_lda, ml_AT, b_ldb, &
                       0.d0, lhs_ml, b_ldc)
            b_ldc = to_blas_int(MSIZE_TDOFS_SCAL)
            b_ldb = to_blas_int(MSIZE_CELL_MAT)
            b_lda = to_blas_int(MSIZE_CELL_SCAL)
            b_m = to_blas_int(gv_cbs)
            b_n = to_blas_int(mk_total_dofs)
            b_k = to_blas_int(mk_gbs_tot)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       1.d0, vm_AT, b_lda, hhoMecaState%grad, b_ldb, &
                       0.d0, lhs_vm, b_ldc)
            b_ldc = to_blas_int(MSIZE_CELL_SCAL)
            b_ldb = to_blas_int(MSIZE_CELL_MAT)
            b_lda = to_blas_int(MSIZE_CELL_SCAL)
            b_m = to_blas_int(gv_cbs)
            b_n = to_blas_int(mk_total_dofs)
            b_k = to_blas_int(mk_gbs_tot)
            call dgemm('N', 'N', b_m, b_n, b_k, &
                       1.d0, lm_AT, b_lda, hhoMecaState%grad, b_ldb, &
                       0.d0, lhs_lm, b_ldc)
! ----- Add stabilization
! ----- += coeff * stab_mk
            do j = 1, mk_total_dofs
                b_n = to_blas_int(mk_total_dofs)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, mk_stab, hhoMecaState%stab(1, j), b_incx, lhs_mm(1, j), &
                           b_incy)
            end do
! ----- += coeff * stab_vv
            do j = 1, gv_total_dofs
                b_n = to_blas_int(gv_total_dofs)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, gv_stab, hhoGVState%stab(1, j), b_incx, lhs_vv(1, j), &
                           b_incy)
            end do
!
! ----- the symmetry is checked inside gdef_log
            lhs_lv = transpose(lhs_vl)
!
! ----- assembly
            call hhoAssGVLhs(hhoCell, hhoData, mapMeca, mapVari, mapLagv, &
                             lhs_mm, lhs_mv, lhs_ml, lhs_vm, lhs_vv, &
                             lhs_vl, lhs_lm, lhs_lv, lhs_ll, lhs)
        end if
!
999     continue
!
! - SYNTHESE DES CODES RETOURS
!
        call codere(cod, hhoQuadCellRigi%nbQuadPoints, hhoComporState%codret)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine check_behavior(hhoComporState)
!
        implicit none
!
        type(HHO_Compor_State), intent(in) :: hhoComporState
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Check behavior
!   In behavior     : type of behavior
! --------------------------------------------------------------------------------------------------
!
        select case (hhoComporState%compor(DEFO))
        case ('GDEF_LOG')
            ASSERT(hhoComporState%l_largestrain)
        case ('PETIT')
            ASSERT(ASTER_FALSE)
            ASSERT(.not. hhoComporState%l_largestrain)
        case default
            ASSERT(ASTER_FALSE)
        end select
        ASSERT(.not. hhoComporState%c_plan)
        ASSERT(.not. hhoComporState%axis)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine gdef_log(BEHinteg, hhoCS, ndim, ipg, time_prev, &
                        time_curr, F_prev, F_curr, var_prev, var_curr, &
                        lag_prev, lag_curr, GV_prev, GV_curr, PK1_curr, &
                        sig_vari, sig_lagv, sig_gv, dPK1_dF, dPK1_dv, &
                        dPK1_dl, dsv_dF, dsv_dv, dsv_dl, dsl_dF, &
                        dsl_dl, dsgv_dgv, cod)
!
        implicit none
!
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        type(HHO_Compor_State), intent(inout) :: hhoCS
        integer, intent(in) :: ndim
        integer, intent(in) :: ipg
        real(kind=8), intent(in) :: time_prev
        real(kind=8), intent(in) :: time_curr
        real(kind=8), intent(in) :: F_prev(3, 3)
        real(kind=8), intent(in) :: F_curr(3, 3)
        real(kind=8), intent(in) :: var_prev, var_curr
        real(kind=8), intent(in) :: lag_prev, lag_curr
        real(kind=8), intent(in) :: GV_prev(3), GV_curr(3)
        real(kind=8), intent(out) :: PK1_curr(3, 3)
        real(kind=8), intent(out) :: sig_vari, sig_lagv, sig_gv(3)
        real(kind=8), intent(out) :: dsv_dv, dsv_dl, dsl_dl, dsgv_dgv(3, 3)
        real(kind=8), intent(out) :: dPK1_dF(3, 3, 3, 3), dPK1_dv(3, 3), dPK1_dl(3, 3)
        real(kind=8), intent(out) :: dsv_dF(3, 3), dsl_dF(3, 3)
        integer, intent(out) :: cod
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the behavior laws for GDEF_LOF
!   IO BEHinteg     : integration informations
!   In ndim         : dimension of the problem
!   In fami         : familly of quadrature points
!   In typmod       : type of modelization
!   In imate        : materiau code
!   In compor       : type of behavior
!   In option       : option of computations
!   In carcri       : local criterion of convergence
!   In lgpg         : size of internal variables for 1 pg
!   In ipg          : i-th quadrature point
!   In time_prev    : previous time T-
!   In time_curr    : current time T+
!   In angmas       : LES TROIS ANGLES DU MOT_CLEF MASSIF
!   In multcomp     : ?
!   In cplan        : plane stress hypothesis
!   In F_prev       : previous deformation gradient at T-
!   In F_curr       : curr deformation gradient at T+
!   In sig_prev_pg  : cauchy stress at T-  (XX, YY, ZZ, XY, XZ, YZ)
!   In vi_prev_pg   : internal variables at T-
!   Out sig_curr    : cauchy stress at T+  (XX, YY, ZZ, XY, XZ, YZ)
!   Out vi_curr     : internal variables at T+
!   Out Pk1_curr    : piola-kirschooff 1 at T+
!   Out module_tang : tangent modulus dPK1/dF(Fp)
!   Out cod         : info on integration of the LDC
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: gn(3, 3), lamb(3), logl(3), epslPrev(6), epslIncr(6)
        real(kind=8) :: tlogPrev(6), tlogCurr(6)
        real(kind=8) :: dtde(6, 6), PK2_prev(6), PK2_curr(6), sig(6)
        real(kind=8) :: dpk2dc(6, 6), ftot(3, 3, 3, 3)
        real(Kind=8) :: eplcm(3*ndim+2), eplci(3*ndim+2)
        real(Kind=8) :: silcm(3*ndim+2), silcp(3*ndim+2), dsde(3*ndim+2, 3*ndim+2)
        real(kind=8) :: sigPrev(11), viPrev(hhoCS%lgpg), viCurr(hhoCS%lgpg)
        real(kind=8) :: dT_dv(6), dT_dl(6)
        real(kind=8) :: dsv_de(6), dsl_de(6)
        real(kind=8) :: dsl_dv, norm, pe(3, 3, 3, 3)
        integer :: lgpg, imate, neu, neg, ntot
        aster_logical :: lCorr, lMatr, lSigm, lVari
!
        lCorr = L_CORR(hhoCS%option)
        lMatr = L_MATR(hhoCS%option)
        lSigm = L_SIGM(hhoCS%option)
        lVari = L_VARI(hhoCS%option)
!
        neu = 2*ndim
        neg = 2+ndim
        ntot = neu+neg
!
        lgpg = hhoCS%lgpg
        imate = hhoCS%imater
        sigPrev(1:ntot) = hhoCS%sig_prev((ipg-1)*ntot+1:ipg*ntot)
        viPrev = hhoCS%vari_prev((ipg-1)*lgpg+1:ipg*lgpg)
!
! ----- Compute pre-processing Elog
!
        call prelog(ndim, lgpg, viPrev, gn, lamb, &
                    logl, F_prev, F_curr, epslPrev, epslIncr, &
                    tlogPrev, lCorr, cod)
        if (cod .ne. 0) then
            goto 999
        end if
! Preparation des deformations generalisees de ldc en t- et t+
        eplcm(1:neu) = epslPrev(1:neu)
        eplci(1:neu) = epslIncr(1:neu)
        eplcm(neu+1) = var_prev
        eplci(neu+1) = var_curr-var_prev
        eplcm(neu+2) = lag_prev
        eplci(neu+2) = lag_curr-lag_prev
        eplcm(neu+2+1:neu+2+ndim) = GV_prev(1:ndim)
        eplci(neu+2+1:neu+2+ndim) = GV_curr(1:ndim)-GV_prev(1:ndim)
! Preparation des contraintes generalisees de ldc en t-
        silcm(1:neu) = viPrev(lgpg-5:lgpg-6+neu)
        silcm(neu+1:ntot) = sigPrev(neu+1:ntot)
!
! ----- Compute Stress and module_tangent
!
        silcp = 0.d0
        viCurr = 0.d0
        dsde = 0.d0
        call nmcomp(BEHinteg, hhoCS%fami, ipg, 1, ndim, &
                    hhoCS%typmod, imate, hhoCS%compor, hhoCS%carcri, time_prev, &
                    time_curr, ntot, eplcm, eplci, ntot, &
                    silcm, viPrev, hhoCS%option, hhoCS%angl_naut, silcp, &
                    viCurr, ntot*ntot, dsde, cod, hhoCS%mult_comp)
!
! ----- Test the code of the LDC
!
        if (cod .eq. 1) goto 999
!
! Archivage des contraintes mecaniques en t+ (tau tilda) dans les vi
        if (lVari) then
            viCurr(lgpg-1:lgpg) = 0.d0
            viCurr(lgpg-5:lgpg-6+neu) = silcp(1:neu)
            hhoCS%vari_curr((ipg-1)*lgpg+1:ipg*lgpg) = viCurr
        end if
!
! ----- Compute post-processing Elog
!
        dtde = 0.d0
        dtde(1:neu, 1:neu) = dsde(1:neu, 1:neu)
        tlogCurr = 0.d0
        tlogCurr(1:neu) = silcp(1:neu)
!
        call poslog(lCorr, lMatr, lSigm, lVari, tlogPrev, &
                    tlogCurr, F_prev, lgpg, viCurr, ndim, &
                    F_curr, ipg, dtde, sigPrev, hhoCS%c_plan, &
                    hhoCS%fami, imate, time_curr, hhoCS%angl_naut, gn, &
                    lamb, logl, sig, dpk2dc, PK2_prev, &
                    PK2_curr, cod)
!
! ----- Test the code of the LDC
!
        if (cod .ne. 0) goto 999
!
        if (.not. lCorr) then
            PK2_curr = PK2_prev
            tlogCurr = tlogPrev
        end if
!
        if (lSigm) then
            hhoCS%sig_curr((ipg-1)*ntot+1:(ipg-1)*ntot+neu) = sig(1:neu)
            hhoCS%sig_curr((ipg-1)*ntot+neu+1:ipg*ntot) = silcp(neu+1:ntot)
        end if
!
! ----- Compute PK1
!
        call deflg4(gn, lamb, logl, F_curr, pe)
        call prodmt(tlogCurr, pe, Pk1_curr)
        sig_vari = silcp(neu+1)
        sig_lagv = silcp(neu+2)
        sig_gv = 0.d0
        sig_gv(1:ndim) = silcp(neu+2+1:ntot)
!
        dPK1_dF = 0.d0
        dPK1_dv = 0.d0
        dPK1_dl = 0.d0
        dsv_dF = 0.d0
        dsv_dv = 0.d0
        dsv_dl = 0.d0
        dsl_dF = 0.d0
        dsl_dl = 0.d0
        dsgv_dgv = 0.d0
!
        if (lMatr) then
!
! ----- Unpack lagrangian tangent modulus
!
            call desymt46(dpk2dc, ftot)
!
! ----- Compute nominal tangent modulus
!
            call lagmodtonommod(ftot, PK2_curr, F_curr, dPK1_dF)
!
            dT_dv = 0.d0
            dT_dl = 0.d0
            dsv_de = 0.d0
            dsl_de = 0.d0
!
            dsv_dv = dsde(neu+1, neu+1)
            dsv_dl = dsde(neu+1, neu+2)
            dsl_dv = dsde(neu+2, neu+1)
            dsl_dl = dsde(neu+2, neu+2)
            dT_dv(1:neu) = dsde(1:neu, neu+1)
            dT_dl(1:neu) = dsde(1:neu, neu+2)
            dsv_de(1:neu) = dsde(neu+1, 1:neu)
            dsl_de(1:neu) = dsde(neu+2, 1:neu)
            dsgv_dgv(1:ndim, 1:ndim) = dsde(neu+3:ntot, neu+3:ntot)
!
!call hhoPrintMat(dsde)
!
!           dP_da = dT_da * dElog_dF
            call prodmt(dT_dv, pe, dPK1_dv)
            call prodmt(dT_dl, pe, dPK1_dl)
!
!           da_dF = da_dElog * dElog_dF
            call prodmt(dsv_de, pe, dsv_dF)
            call prodmt(dsl_de, pe, dsl_dF)
!
! ----- Verify symmetry
            norm = max(1.d0, dsv_dl)
            ASSERT(abs(dsv_dl-dsl_dv) < 1d-8*norm)
        end if
! print *, hhoCS%option
! print *, ipg, lgpg, ntot
! print *, eplcm
! print *, eplci
! print *, sig_vari, sig_lagv, sig_gv
! print *, dsv_dv, dsv_dl, dsl_dl, dsgv_dgv
! print*, ipg, dsv_dv, dsv_dl, dsl_dl
! print*, dPK1_dv
! print*, dPK1_dl
! print*, dsv_dF
! print*, dsl_dF
!if (abs(dsl_dl) < 1d-8) dsl_dl = 1.d0
!
999     continue
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine petit(BEHinteg, hhoCS, ndim, ipg, time_prev, &
                     time_curr, Eps_prev, Eps_curr, var_prev, var_curr, &
                     lag_prev, lag_curr, GV_prev, GV_curr, Sig_curr, &
                     sig_vari, sig_lagv, sig_gv, dSig_dEps, dSig_dv, &
                     dSig_dl, dsv_dEps, dsv_dv, dsv_dl, dsl_dEps, &
                     dsl_dl, dsgv_dgv, cod)
!
        implicit none
!
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        type(HHO_Compor_State), intent(inout) :: hhoCS
        integer, intent(in) :: ndim
        integer, intent(in) :: ipg
        real(kind=8), intent(in) :: time_prev
        real(kind=8), intent(in) :: time_curr
        real(kind=8), intent(in) :: Eps_prev(6)
        real(kind=8), intent(in) :: Eps_curr(6)
        real(kind=8), intent(in) :: var_prev, var_curr
        real(kind=8), intent(in) :: lag_prev, lag_curr
        real(kind=8), intent(in) :: GV_prev(3), GV_curr(3)
        real(kind=8), intent(out) :: Sig_curr(6)
        real(kind=8), intent(out) :: sig_vari, sig_lagv, sig_gv(3)
        real(kind=8), intent(out) :: dsv_dv, dsv_dl, dsl_dl, dsgv_dgv(3, 3)
        real(kind=8), intent(out) :: dSig_dEps(6, 6), dSig_dv(6), dSig_dl(6)
        real(kind=8), intent(out) :: dsv_dEps(6), dsl_dEps(6)
        integer, intent(out) :: cod
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the behavior laws for GDEF_LOF
!   IO BEHinteg     : integration informations
!   In ndim         : dimension of the problem
!   In fami         : familly of quadrature points
!   In typmod       : type of modelization
!   In imate        : materiau code
!   In compor       : type of behavior
!   In option       : option of computations
!   In carcri       : local criterion of convergence
!   In lgpg         : size of internal variables for 1 pg
!   In ipg          : i-th quadrature point
!   In time_prev    : previous time T-
!   In time_curr    : current time T+
!   In angmas       : LES TROIS ANGLES DU MOT_CLEF MASSIF
!   In multcomp     : ?
!   In cplan        : plane stress hypothesis
!   In F_prev       : previous deformation gradient at T-
!   In F_curr       : curr deformation gradient at T+
!   In sig_prev_pg  : cauchy stress at T-  (XX, YY, ZZ, XY, XZ, YZ)
!   In vi_prev_pg   : internal variables at T-
!   Out sig_curr    : cauchy stress at T+  (XX, YY, ZZ, XY, XZ, YZ)
!   Out vi_curr     : internal variables at T+
!   Out Pk1_curr    : piola-kirschooff 1 at T+
!   Out module_tang : tangent modulus dPK1/dF(Fp)
!   Out cod         : info on integration of the LDC
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: Cauchy_prev(6), Cauchy_curr(6), sig(6)
        real(Kind=8) :: eplcm(3*ndim+2), eplci(3*ndim+2)
        real(Kind=8) :: silcm(3*ndim+2), silcp(3*ndim+2), dsde(3*ndim+2, 3*ndim+2)
        real(kind=8) :: sigPrev(11), viPrev(hhoCS%lgpg), viCurr(hhoCS%lgpg)
        real(kind=8) :: dsl_dv
        integer :: lgpg, imate, neu, neg, ntot
        aster_logical :: lCorr, lMatr, lSigm, lVari
!
        lCorr = L_CORR(hhoCS%option)
        lMatr = L_MATR(hhoCS%option)
        lSigm = L_SIGM(hhoCS%option)
        lVari = L_VARI(hhoCS%option)
!
        neu = 2*ndim
        neg = 2+ndim
        ntot = neu+neg
!
        lgpg = hhoCS%lgpg
        imate = hhoCS%imater
        sigPrev(1:ntot) = hhoCS%sig_prev((ipg-1)*ntot+1:ipg*ntot)
        call tranfoMatToSym(ndim, sigPrev(1:neu), Cauchy_prev)
        viPrev = hhoCS%vari_prev((ipg-1)*lgpg+1:ipg*lgpg)
! Preparation des deformations generalisees de ldc en t- et t+
        eplcm(1:neu) = Eps_prev(1:neu)
        eplci(1:neu) = Eps_curr(1:neu)-Eps_prev(1:neu)
        eplcm(neu+1) = var_prev
        eplci(neu+1) = var_curr-var_prev
        eplcm(neu+2) = lag_prev
        eplci(neu+2) = lag_curr-lag_prev
        eplcm(neu+2+1:neu+2+ndim) = GV_prev(1:ndim)
        eplci(neu+2+1:neu+2+ndim) = GV_curr(1:ndim)-GV_prev(1:ndim)
! Preparation des contraintes generalisees de ldc en t-
        silcm(1:neu) = Cauchy_prev(1:neu)
        silcm(neu+1:ntot) = sigPrev(neu+1:ntot)
!
! ----- Compute Stress and module_tangent
!
        silcp = 0.d0
        viCurr = 0.d0
        dsde = 0.d0
        call nmcomp(BEHinteg, hhoCS%fami, ipg, 1, ndim, &
                    hhoCS%typmod, imate, hhoCS%compor, hhoCS%carcri, time_prev, &
                    time_curr, ntot, eplcm, eplci, ntot, &
                    silcm, viPrev, hhoCS%option, hhoCS%angl_naut, silcp, &
                    viCurr, ntot*ntot, dsde, cod, hhoCS%mult_comp)
!
! ----- Test the code of the LDC
!
        if (cod .eq. 1) goto 999
!
! Archivage des contraintes mecaniques en t+ (tau tilda) dans les vi
        if (lVari) then
            hhoCS%vari_curr((ipg-1)*lgpg+1:ipg*lgpg) = viCurr
        end if
!
! ----- Compute post-processing Elog
!
        Cauchy_curr = 0.d0
        Cauchy_curr(1:neu) = silcp(1:neu)
!
!
! --------- For new prediction and nmisot.F90
        if (L_PRED(hhoCS%option)) then
            Cauchy_curr = 0.d0
        end if
!
        if (.not. lCorr) then
            Cauchy_curr = Cauchy_prev
        end if
!
        if (lSigm) then
            call tranfoSymToMat(ndim, Cauchy_curr, sig)
            hhoCS%sig_curr((ipg-1)*ntot+1:(ipg-1)*ntot+neu) = sig(1:neu)
            hhoCS%sig_curr((ipg-1)*ntot+neu+1:ipg*ntot) = silcp(neu+1:ntot)
        end if
!
! ----- Compute stress
!
        Sig_curr = 0.d0
        Sig_curr(1:neu) = Cauchy_curr(1:neu)
        sig_vari = silcp(neu+1)
        sig_lagv = silcp(neu+2)
        sig_gv = 0.d0
        sig_gv(1:ndim) = silcp(neu+2+1:ntot)
!
        dSig_dEps = 0.d0
        dSig_dv = 0.d0
        dSig_dl = 0.d0
        dsv_dEps = 0.d0
        dsv_dv = 0.d0
        dsv_dl = 0.d0
        dsl_dEps = 0.d0
        dsl_dl = 0.d0
        dsgv_dgv = 0.d0
!
        if (lMatr) then
!
            dSig_dEps(1:neu, 1:neu) = dsde(1:neu, 1:neu)
            dsv_dv = dsde(neu+1, neu+1)
            dsv_dl = dsde(neu+1, neu+2)
            dsl_dv = dsde(neu+2, neu+1)
            dsl_dl = dsde(neu+2, neu+2)
            dSig_dv(1:neu) = dsde(1:neu, neu+1)
            dSig_dl(1:neu) = dsde(1:neu, neu+2)
            dsv_dEps(1:neu) = dsde(neu+1, 1:neu)
            dsl_dEps(1:neu) = dsde(neu+2, 1:neu)
            dsgv_dgv(1:ndim, 1:ndim) = dsde(neu+3:ntot, neu+3:ntot)
!
! print *, hhoCS%option
! print *, ipg, lgpg, ntot
! print *, eplcm
! print *, eplci
! print *, sig_vari, sig_lagv, sig_gv
! print *, dsv_dv, dsv_dl, dsl_dl, dsgv_dgv
!if (abs(dsl_dl) < 1d-8) dsl_dl = 1.d0
        end if
!
999     continue
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine numGVMap(hhoCell, hhoData, mapMeca, mapVari, mapLagv)
!
        implicit none
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        integer, intent(out) :: mapMeca(MSIZE_TDOFS_VEC)
        integer, intent(out) :: mapVari(MSIZE_TDOFS_SCAL)
        integer, intent(out) :: mapLagv(MSIZE_CELL_SCAL)
!
!--------------------------------------------------------------------------------------------------
! Numbering map for GRAD_VARI
!
!--------------------------------------------------------------------------------------------------
        integer :: mk_cbs, mk_fbs, mk_total_dofs, gv_cbs, gv_fbs, gv_total_dofs
        integer :: i_face, i_dof, num_tot, num_gv, num_vari, num_mk
!
!
        call hhoMecaDofs(hhoCell, hhoData, mk_cbs, mk_fbs, mk_total_dofs)
        call hhoTherDofs(hhoCell, hhoData, gv_cbs, gv_fbs, gv_total_dofs)
!
        mapMeca = 0
        mapVari = 0
        mapLagv = 0
!
        num_tot = 0
        num_gv = 0
        num_mk = 0
        num_vari = 0
!
        do i_face = 1, hhoCell%nbfaces
            do i_dof = 1, mk_fbs
                num_tot = num_tot+1
                num_mk = num_mk+1
                mapMeca(mk_cbs+num_mk) = num_tot
            end do
            do i_dof = 1, gv_fbs
                num_tot = num_tot+1
                num_gv = num_gv+1
                mapVari(gv_cbs+num_gv) = num_tot
            end do
        end do
!
        do i_dof = 1, mk_cbs
            num_tot = num_tot+1
            mapMeca(i_dof) = num_tot
        end do
        do i_dof = 1, gv_cbs
            num_tot = num_tot+1
            mapVari(i_dof) = num_tot
        end do
        do i_dof = 1, gv_cbs
            num_tot = num_tot+1
            mapLagv(i_dof) = num_tot
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoAssGVRhs(hhoCell, hhoData, mapMeca, mapVari, mapLagv, &
                           rhs_meca, rhs_vari, rhs_lagv, rhs)
!
        implicit none
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        integer, intent(in) :: mapMeca(MSIZE_TDOFS_VEC)
        integer, intent(in) :: mapVari(MSIZE_TDOFS_SCAL)
        integer, intent(in) :: mapLagv(MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: rhs_meca(MSIZE_TDOFS_VEC)
        real(kind=8), intent(in) :: rhs_vari(MSIZE_TDOFS_SCAL)
        real(kind=8), intent(in) :: rhs_lagv(MSIZE_CELL_SCAL)
        real(kind=8), intent(out) :: rhs(MSIZE_TDOFS_MIX)
!--------------------------------------------------------------------------------------------------
! Assembly RHS for GRAD_VARI
!
!--------------------------------------------------------------------------------------------------
        integer :: mk_cbs, mk_fbs, mk_total_dofs, gv_cbs, gv_fbs, gv_total_dofs
        integer :: i_dof
!
        call hhoMecaDofs(hhoCell, hhoData, mk_cbs, mk_fbs, mk_total_dofs)
        call hhoTherDofs(hhoCell, hhoData, gv_cbs, gv_fbs, gv_total_dofs)
!
        rhs = 0.d0
!
        do i_dof = 1, mk_total_dofs
            rhs(mapMeca(i_dof)) = rhs_meca(i_dof)
        end do
        do i_dof = 1, gv_total_dofs
            rhs(mapVari(i_dof)) = rhs_vari(i_dof)
        end do
        do i_dof = 1, gv_cbs
            rhs(mapLagv(i_dof)) = rhs_lagv(i_dof)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoAssGVLhs(hhoCell, hhoData, mapMeca, mapVari, mapLagv, &
                           lhs_mm, lhs_mv, lhs_ml, lhs_vm, lhs_vv, &
                           lhs_vl, lhs_lm, lhs_lv, lhs_ll, lhs)
!
        implicit none
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        integer, intent(in) :: mapMeca(MSIZE_TDOFS_VEC)
        integer, intent(in) :: mapVari(MSIZE_TDOFS_SCAL)
        integer, intent(in) :: mapLagv(MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: lhs_mm(MSIZE_TDOFS_VEC, MSIZE_TDOFS_VEC)
        real(kind=8), intent(in) :: lhs_mv(MSIZE_TDOFS_VEC, MSIZE_TDOFS_SCAL)
        real(kind=8), intent(in) :: lhs_ml(MSIZE_TDOFS_VEC, MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: lhs_vm(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_VEC)
        real(kind=8), intent(in) :: lhs_vv(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL)
        real(kind=8), intent(in) :: lhs_vl(MSIZE_TDOFS_SCAL, MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: lhs_lm(MSIZE_CELL_SCAL, MSIZE_TDOFS_VEC)
        real(kind=8), intent(in) :: lhs_lv(MSIZE_CELL_SCAL, MSIZE_TDOFS_SCAL)
        real(kind=8), intent(in) :: lhs_ll(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
        real(kind=8), intent(out) :: lhs(MSIZE_TDOFS_MIX, MSIZE_TDOFS_MIX)
!--------------------------------------------------------------------------------------------------
! Assembly LHS for GRAD_VARI
!
!--------------------------------------------------------------------------------------------------
        integer :: mk_cbs, mk_fbs, mk_total_dofs, gv_cbs, gv_fbs, gv_total_dofs
        integer :: i_row, i_col
!
        call hhoMecaDofs(hhoCell, hhoData, mk_cbs, mk_fbs, mk_total_dofs)
        call hhoTherDofs(hhoCell, hhoData, gv_cbs, gv_fbs, gv_total_dofs)
!
        lhs = 0.d0
!
! --- Bloc meca
        do i_row = 1, mk_total_dofs
            do i_col = 1, mk_total_dofs
                lhs(mapMeca(i_row), mapMeca(i_col)) = lhs_mm(i_row, i_col)
            end do
            do i_col = 1, gv_total_dofs
                lhs(mapMeca(i_row), mapVari(i_col)) = lhs_mv(i_row, i_col)
            end do
            do i_col = 1, gv_cbs
                lhs(mapMeca(i_row), mapLagv(i_col)) = lhs_ml(i_row, i_col)
            end do
        end do
!
! --- Bloc vari
        do i_row = 1, gv_total_dofs
            do i_col = 1, mk_total_dofs
                lhs(mapVari(i_row), mapMeca(i_col)) = lhs_vm(i_row, i_col)
            end do
            do i_col = 1, gv_total_dofs
                lhs(mapVari(i_row), mapVari(i_col)) = lhs_vv(i_row, i_col)
            end do
            do i_col = 1, gv_cbs
                lhs(mapVari(i_row), mapLagv(i_col)) = lhs_vl(i_row, i_col)
            end do
        end do
!
! --- Bloc lagr
        do i_row = 1, gv_cbs
            do i_col = 1, mk_total_dofs
                lhs(mapLagv(i_row), mapMeca(i_col)) = lhs_lm(i_row, i_col)
            end do
            do i_col = 1, gv_total_dofs
                lhs(mapLagv(i_row), mapVari(i_col)) = lhs_lv(i_row, i_col)
            end do
            do i_col = 1, gv_cbs
                lhs(mapLagv(i_row), mapLagv(i_col)) = lhs_ll(i_row, i_col)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_gv(this, hhoCell, hhoData, hhoComporState)
!
        implicit none
!
        class(HHO_GV_State), intent(inout) :: this
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(inout) :: hhoData
        type(HHO_Compor_State), intent(in) :: hhoComporState
!
! --------------------------------------------------------------------------------------------------
!
!  initialize HHO_GV_STATE
! --------------------------------------------------------------------------------------------------
!
        integer :: num_tot, num_gv, iFace, idof
        integer :: mk_cbs, mk_fbs, mk_total_dofs
        integer :: gv_cbs, gv_fbs, gv_total_dofs, total_dofs
        real(kind=8) :: tmp_prev(MSIZE_TDOFS_MIX), tmp_incr(MSIZE_TDOFS_MIX)
        aster_logical :: forc_noda
        blas_int :: b_incx, b_incy, b_n
!
        forc_noda = hhoComporState%option == "FORC_NODA"
        if (hhoComporState%option .ne. "RIGI_MECA") then
            if (hhoComporState%typmod(2) .ne. "GRADVARI") then
                ASSERT(ASTER_FALSE)
            end if
!
            call hhoMecaDofs(hhoCell, hhoData, mk_cbs, mk_fbs, mk_total_dofs)
            call hhoTherDofs(hhoCell, hhoData, gv_cbs, gv_fbs, gv_total_dofs)
            total_dofs = mk_total_dofs+gv_total_dofs+gv_cbs
!
            call readVector('PDEPLMR', total_dofs, tmp_prev)
            if (.not. forc_noda) then
                call readVector('PDEPLPR', total_dofs, tmp_incr)
            end if
!
            num_tot = 0
            num_gv = 0
            do iFace = 1, hhoCell%nbfaces
                num_tot = num_tot+mk_fbs
                do iDof = 1, gv_fbs
                    num_tot = num_tot+1
                    num_gv = num_gv+1
                    this%vari_prev(gv_cbs+num_gv) = tmp_prev(num_tot)
                    if (.not. forc_noda) then
                        this%vari_incr(gv_cbs+num_gv) = tmp_incr(num_tot)
                    end if
                end do
            end do
            num_tot = num_tot+mk_cbs
            do iDof = 1, gv_cbs
                num_tot = num_tot+1
                this%vari_prev(iDof) = tmp_prev(num_tot)
                if (.not. forc_noda) then
                    this%vari_incr(iDof) = tmp_incr(num_tot)
                end if
            end do
            do iDof = 1, gv_cbs
                num_tot = num_tot+1
                this%lagv_prev(iDof) = tmp_prev(num_tot)
                if (.not. forc_noda) then
                    this%lagv_incr(iDof) = tmp_incr(num_tot)
                end if
            end do
        end if
!
! --- compute in T+
!
        b_n = to_blas_int(gv_total_dofs)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, this%vari_prev, b_incx, this%vari_curr, b_incy)
        b_n = to_blas_int(gv_cbs)
        call dcopy(b_n, this%lagv_prev, b_incx, this%lagv_curr, b_incy)
!
        if (.not. forc_noda) then
            b_n = to_blas_int(gv_total_dofs)
            call daxpy(b_n, 1.d0, this%vari_incr, b_incx, this%vari_curr, &
                       b_incy)
            b_n = to_blas_int(gv_cbs)
            call daxpy(b_n, 1.d0, this%lagv_incr, b_incx, this%lagv_curr, &
                       b_incy)

        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalcOpGv(hhoCell, hhoData, l_largestrains, hhoMecaState, hhoGvState)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        aster_logical, intent(in) :: l_largestrains
        type(HHO_Meca_State), intent(inout) :: hhoMecaState
        type(HHO_GV_State), intent(inout) :: hhoGvState
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Mechanic
!
! Compute operators for mechanic
!
! --------------------------------------------------------------------------------------------------
!
! In  hhoCell         : hho Cell
! In hhoData          : information about the HHO formulation
! In l_largestrains   : large strains ?
! Out gradfull        : full gradient for mechanics
! Out stab            : stabilization for mechanics
! --------------------------------------------------------------------------------------------------
!
        integer :: jgrad, jstab, gv_cbs, gv_fbs, gv_total_dofs, gv_gbs, j
        blas_int :: b_incx, b_incy, b_n
!
        if (ASTER_FALSE) then
            call jevech('PCHHOGT', 'L', jgrad)
            call jevech('PCHHOST', 'L', jstab)
!
            call hhoReloadPreCalcMeca(hhoCell, hhoData, l_largestrains, zr(jgrad), zr(jstab), &
                                      hhoMecaState%grad, hhoMecaState%stab)
!
            call hhoTherNLDofs(hhoCell, hhoData, gv_cbs, gv_fbs, gv_total_dofs, &
                               gv_gbs)
            if (l_largestrains) then
                do j = 1, gv_total_dofs
                    b_n = to_blas_int(gv_gbs)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(jgrad+(j-1)*gv_gbs), b_incx, hhoGVState%grad(1, j), &
                               b_incy)
                end do
            else
                call hhoCalcOpTher(hhoCell, hhoData, hhoGVState%grad)
            end if
            do j = 1, gv_total_dofs
                b_n = to_blas_int(gv_total_dofs)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(jstab+(j-1)*gv_total_dofs), b_incx, hhoGVState%stab(1, j), &
                           b_incy)
            end do
        else
            call hhoCalcOpMeca(hhoCell, hhoData, l_largestrains, hhoMecaState%grad, &
                               hhoMecaState%stab)
            call hhoCalcOpTher(hhoCell, hhoData, hhoGVState%grad, hhoGVState%stab)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsLargeMG(hhoCell, module_tang, weight, BSCEval, gv_cbs, &
                                    BSCEvalG, mk_gbs, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: module_tang(3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: BSCEvalG(MSIZE_CELL_SCAL)
        integer, intent(in) :: gv_cbs, mk_gbs
        real(kind=8), intent(inout) :: AT(MSIZE_CELL_MAT, MSIZE_CELL_SCAL)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (module_tang:cphi, gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Acphi(3, 3, MSIZE_CELL_SCAL)
        integer :: i, j, k, row, gbs_cmp
        blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------------------------------------------
!
! --------- Eval (module_tang : scphi)_T
        do i = 1, gv_cbs
            qp_Acphi(:, :, i) = weight*module_tang*BSCEval(i)
        end do
!
! -------- Compute scalar_product of (C_sgphi(j), gphi(j))_T
        gbs_cmp = mk_gbs/(hhoCell%ndim*hhoCell%ndim)
!
        row = 1
        do i = 1, hhoCell%ndim
            do j = 1, hhoCell%ndim
                do k = 1, gbs_cmp
                    b_n = to_blas_int(gv_cbs)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, BSCEvalG(k), qp_Acphi(i, j, :), b_incx, AT(row, :), &
                               b_incy)
                    row = row+1
                end do
            end do
        end do
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsLargeGM(hhoCell, module_tang, weight, BSCEval, gv_cbs, &
                                    BSCEvalG, mk_gbs, AT)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        real(kind=8), intent(in) :: module_tang(3, 3)
        real(kind=8), intent(in) :: weight
        real(kind=8), intent(in) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8), intent(in) :: BSCEvalG(MSIZE_CELL_SCAL)
        integer, intent(in) :: gv_cbs, mk_gbs
        real(kind=8), intent(inout) :: AT(MSIZE_CELL_SCAL, MSIZE_CELL_MAT)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product AT += (cphi, module_tang:gphi)_T at a quadrature point
!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto-plastic tangent moduli
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In gbs_cmp      : size of BSCEval
!   In gbs          : number of rows of AT
!   Out AT          : contribution of At
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_Agphi(3, 3, MSIZE_CELL_SCAL)
        integer :: i, j, k, col, gbs_cmp
        blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------------------------------------------
!
! --------- Eval (module_tang : sgphi)_T
        gbs_cmp = mk_gbs/(hhoCell%ndim*hhoCell%ndim)
!
        do i = 1, gbs_cmp
            qp_Agphi(:, :, i) = weight*module_tang*BSCEvalG(i)
        end do
!
! -------- Compute scalar_product of (C_sgphi(j), gphi(j))_T
        col = 1
        do i = 1, hhoCell%ndim
            do j = 1, hhoCell%ndim
                do k = 1, gbs_cmp
                    b_n = to_blas_int(gv_cbs)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, qp_Agphi(i, j, k), BSCEval, b_incx, AT(:, col), &
                               b_incy)
                    col = col+1
                end do
            end do
        end do
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine forc_noda_stress(hhoCS, ndim, ipg, F_curr, PK1_curr, &
                                Cauchy_curr, sig_vari, sig_lagv, sig_gv)
!
        implicit none
!
        type(HHO_Compor_State), intent(inout) :: hhoCS
        integer, intent(in) :: ndim
        integer, intent(in) :: ipg
        real(kind=8), intent(in) :: F_curr(3, 3)
        real(kind=8), intent(out) :: PK1_curr(3, 3), Cauchy_curr(6)
        real(kind=8), intent(out) :: sig_vari, sig_lagv, sig_gv(3)
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Get stress for FORC_NODA
! --------------------------------------------------------------------------------------------------
!
        integer :: neu, neg, ntot
        real(kind=8) :: sigPrev(11)
! --------------------------------------------------------------------------------------------------
!
        neu = 2*ndim
        neg = 2+ndim
        ntot = neu+neg
!
        sigPrev(1:ntot) = hhoCS%sig_prev((ipg-1)*ntot+1:ipg*ntot)
        sig_vari = sigPrev(neu+1)
        sig_lagv = sigPrev(neu+2)
        sig_gv = 0.d0
        sig_gv(1:ndim) = sigPrev(neu+2+1:ntot)
!
        call tranfoMatToSym(ndim, sigPrev(1:neu), Cauchy_curr)
!
        PK1_curr = 0.d0
        if (hhoCS%l_largestrain) then
            call sigtopk1(ndim, Cauchy_curr, F_curr, PK1_curr)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    real(kind=8) function hhoCalcStabCoeffGV(fami, npg)
!
        implicit none
!
        character(len=4) :: fami
        integer, intent(in) :: npg
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  GRAD_VARI - Evaluate stabilzation coefficient
! --------------------------------------------------------------------------------------------------
!
! --- Local variables
!
        integer :: jmate, imate
        integer :: ipg, iok(1)
        real(kind=8) :: vale(1)
!
        call jevech('PMATERC', 'L', jmate)
        imate = zi(jmate-1+1)
        hhoCalcStabCoeffGV = 0.d0
!
        do ipg = 1, npg
            call rcvalb(fami, ipg, 1, '+', imate, &
                        ' ', 'NON_LOCAL', 0, ' ', [0.d0], &
                        1, ['C_GRAD_VARI'], vale, iok, 1)
            hhoCalcStabCoeffGV = hhoCalcStabCoeffGV+vale(1)
        end do
!
        hhoCalcStabCoeffGV = 10.d0*hhoCalcStabCoeffGV/real(npg, kind=8)
        ASSERT(hhoCalcStabCoeffGV > 0.d0)
!
    end function
!
end module
