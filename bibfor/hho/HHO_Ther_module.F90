! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
module HHO_Ther_module
!
use HHO_basis_module
use HHO_Dirichlet_module
use HHO_eval_module
use HHO_gradrec_module, only : hhoGradRecFullVec, hhoGradRecVec
use HHO_quadrature_module
use HHO_size_module
use HHO_stabilization_module, only : hhoStabScal, hdgStabScal
use HHO_type
use HHO_utils_module
use NonLin_Datastructure_type
!
implicit none
!
private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/readVector.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
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
    public :: hhoLocalRigiTher, hhoCalcStabCoeff, hhoLocalMassTher
    public :: hhoCalcOpTher
    private :: hhoComputeRhsRigiTher, hhoComputeLhsRigiTher, hhoComputeAgphi
    private :: hhoComputeBehaviourTher, LambdaMax, hhoComputeRhoCpTher
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLocalRigiTher(hhoCell, hhoData, hhoQuadCellRigi, gradrec, stab, &
                                  fami, lhs, rhs)
!
    implicit none
!
        type(HHO_Cell), intent(in)      :: hhoCell
        type(HHO_Data), intent(inout)   :: hhoData
        type(HHO_Quadrature), intent(in):: hhoQuadCellRigi
        real(kind=8), intent(in)        :: gradrec(MSIZE_CELL_VEC, MSIZE_TDOFS_SCAL)
        real(kind=8), intent(in)        :: stab(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL)
        character(len=*), intent(in)    :: fami
        real(kind=8), intent(out), optional :: lhs(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL)
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
        integer:: cbs, fbs, total_dofs, j, faces_dofs, gbs
        integer :: jmate, ipg, icodre(3), jtemps
        real(kind=8), dimension(MSIZE_CELL_VEC) :: bT, G_curr_coeff
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: temp_curr
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8) :: AT(MSIZE_CELL_VEC, MSIZE_CELL_VEC)
        real(kind=8) :: TMP(MSIZE_CELL_VEC, MSIZE_TDOFS_VEC)
        real(kind=8) :: module_tang(3,3), G_curr(3), sig_curr(3)
        real(kind=8) :: coorpg(3), weight, time_curr
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
        call jevech('PTEMPSR', 'L', jtemps)
        time_curr = zr(jtemps)
!
! --- number of dofs
!
        call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)
        faces_dofs = total_dofs - cbs
!
! -- initialization
!
        if(l_lhs) lhs = 0.d0
        if(l_rhs) rhs = 0.d0
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
        if(l_rhs) then
            call readVector('PTEMPER', total_dofs, temp_curr)
            call hhoRenumTherVecInv(hhoCell, hhoData, temp_curr)
        end if
        poum = "+"
!
! ----- compute G_curr = gradrec * temp_curr
!
        call dgemv('N', gbs, total_dofs, 1.d0, gradrec, MSIZE_CELL_VEC, temp_curr, 1,&
                    0.d0, G_curr_coeff,1)
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCellRigi%nbQuadPoints
            coorpg(1:3) = hhoQuadCellRigi%points(1:3,ipg)
            weight = hhoQuadCellRigi%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
            call hhoBasisCell%BSEval(hhoCell, coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
! --------- Eval gradient at T+
!
            G_curr = hhoEvalVecCell(hhoCell, hhoBasisCell, hhoData%grad_degree(), &
                                    coorpg(1:3), G_curr_coeff, gbs)
!
! ------- Compute behavior
!
           call hhoComputeBehaviourTher(phenom, fami, poum, ipg, hhoCell%ndim, time_curr, jmate, &
                                        G_curr, sig_curr, module_tang)
!
! ------- Compute rhs
!
            if(l_rhs) call hhoComputeRhsRigiTher(hhoCell, sig_curr, weight, BSCEval, gbs, bT)
!
! ------- Compute lhs
!
            if(l_lhs) call hhoComputeLhsRigiTher(hhoCell, module_tang, weight, BSCEval, gbs, AT)
!
        end do
!
! ----- compute rhs += Gradrec**T * bT
!
        if(l_rhs) then
            call dgemv('T', gbs, total_dofs, 1.d0, gradrec, MSIZE_CELL_VEC, &
                        bT, 1, 1.d0, rhs, 1)
        end if
!
! ----- compute lhs += gradrec**T * AT * gradrec
! ----- step1: TMP = AT * gradrec
!
        if(l_lhs) then
            call dgemm('N', 'N', gbs, total_dofs, total_dofs, 1.d0, AT, MSIZE_CELL_VEC, &
                   gradrec, MSIZE_CELL_VEC, 0.d0, TMP, MSIZE_CELL_VEC)
!
! ----- step2: lhs += gradrec**T * TMP
!
            call dgemm('T', 'N', total_dofs, total_dofs, gbs, 1.d0, gradrec, MSIZE_CELL_VEC, &
                   TMP, MSIZE_CELL_VEC, 1.d0, lhs, MSIZE_TDOFS_SCAL)
!
        end if
!
! --- add stabilization
!
        call hhoCalcStabCoeff(hhoData, fami, hhoQuadCellRigi%nbQuadPoints, hhoCell%ndim, time_curr)
!
        if(l_rhs) then
                call dsymv('U', total_dofs, hhoData%coeff_stab(), stab, MSIZE_TDOFS_SCAL,&
                   temp_curr, 1, 1.d0, rhs,1)
        end if
!
        if(l_lhs) then
            do j = 1, total_dofs
                call daxpy(total_dofs, hhoData%coeff_stab(), stab(1,j), 1, lhs(1,j), 1)
            end do
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoLocalMassTher(hhoCell, hhoData, hhoQuadCellMass, &
                                  fami,  lhs, rhs)
!
    implicit none
!
        type(HHO_Cell), intent(in)      :: hhoCell
        type(HHO_Data), intent(inout)   :: hhoData
        type(HHO_Quadrature), intent(in):: hhoQuadCellMass
        character(len=*), intent(in)    :: fami
        real(kind=8), intent(out), optional :: lhs(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL)
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
        integer:: cbs, fbs, total_dofs, faces_dofs, gbs
        integer :: jmate, ipg, icodre(3), jtemps
        real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: temp_T_curr
        real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8) :: coorpg(3), weight, time_curr, cp, temp_eval
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
        call jevech('PTEMPSR', 'L', jtemps)
        time_curr = zr(jtemps)
!
! --- number of dofs
!
        call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)
        faces_dofs = total_dofs - cbs
!
! -- initialization
!
        if(l_lhs) lhs = 0.d0
        if(l_rhs) rhs = 0.d0
!
        call hhoBasisCell%initialize(hhoCell)
!
! --- compute temp in T+
!
        temp_T_curr = 0.d0
        if(l_rhs) call readVector('PTEMPER', cbs, temp_T_curr, total_dofs-cbs)
        poum = "+"
!
! ----- Loop on quadrature point
!
        do ipg = 1, hhoQuadCellMass%nbQuadPoints
            coorpg(1:3) = hhoQuadCellMass%points(1:3,ipg)
            weight = hhoQuadCellMass%weights(ipg)
!
! --------- Eval basis function at the quadrature point
!
            call hhoBasisCell%BSEval(hhoCell, coorpg(1:3), 0, hhoData%cell_degree(), BSCEval)
!
! --------- Eval gradient at T+
!
            temp_eval = hhoEvalScalCell(hhoCell, hhoBasisCell, hhoData%cell_degree(),&
                        coorpg, temp_T_curr, cbs)
!
! -------- Compute behavior
!
            call hhoComputeRhoCpTher(phenom, fami, poum, ipg, time_curr, jmate, cp)
!
! -------- Compute rhs
!
            if(l_rhs) call hhoComputeRhsMassTher(temp_eval, cp, weight, BSCEval, cbs, rhs)
!
! -------- Compute lhs
!
            if(l_lhs) call hhoComputeLhsMassTher(cp, weight, BSCEval, cbs, lhs)
!
        end do
!
! ----- Copy the lower part
!
        if(l_lhs )call hhoCopySymPartMat('U', lhs(1:cbs, 1:cbs))
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
        real(kind=8), dimension(MSIZE_CELL_VEC, MSIZE_TDOFS_SCAL), intent(out)   :: gradfull
        real(kind=8), dimension(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL), intent(out), optional  :: stab
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
        real(kind=8), dimension(MSIZE_CELL_SCAL, MSIZE_TDOFS_SCAL) :: gradrec_scal
!
! ----- Compute Gradient reconstruction
        call hhoGradRecFullVec(hhoCell, hhoData, gradfull)
!
! ----- Compute Stabilizatiion
        if(present(stab)) then
            if (hhoData%cell_degree() <= hhoData%face_degree()) then
                call hhoGradRecVec(hhoCell, hhoData, gradrec_scal)
                call hhoStabScal(hhoCell, hhoData, gradrec_scal, stab)
            else if (hhoData%cell_degree() == (hhoData%face_degree() + 1)) then
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
        type(HHO_Cell), intent(in)      :: hhoCell
        real(kind=8), intent(in)        :: stress(3)
        real(kind=8), intent(in)        :: weight
        real(kind=8), intent(in)        :: BSCEval(MSIZE_CELL_SCAL)
        integer, intent(in)             :: gbs
        real(kind=8), intent(inout)     :: bT(MSIZE_CELL_VEC)
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
        real(kind = 8) :: qp_stress(3)
        integer:: i, gbs_cmp, deca
! --------------------------------------------------------------------------------------------------
!
        gbs_cmp = gbs / hhoCell%ndim
        qp_stress = weight * stress
! -------- Compute scalar_product of (stress, gphi)_T
! -------- (RAPPEL: the composents of the gradient are saved by rows)
        deca = 0
        do i = 1, hhoCell%ndim
            call daxpy(gbs_cmp, qp_stress(i), BSCEval, 1, bT(deca+1), 1)
            deca = deca + gbs_cmp
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeRhsMassTher(temp_curr, cp, weight, BSCEval, cbs, rhs)
!
    implicit none
!
        real(kind=8), intent(in)        :: temp_curr, cp
        real(kind=8), intent(in)        :: weight
        real(kind=8), intent(in)        :: BSCEval(MSIZE_CELL_SCAL)
        integer, intent(in)             :: cbs
        real(kind=8), intent(inout)     :: rhs(MSIZE_TDOFS_SCAL)
!
! ------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the scalar product rhs += (rho_cp * temp, vT)_T at a quadrature point
!   In weight       : quadrature weight
!   In BSCEval      : Basis of one composant gphi
!   In cbs          : number of rows of bT
!   Out bT          : contribution of bt
! --------------------------------------------------------------------------------------------------
!
        real(kind = 8) :: qp_temp
! --------------------------------------------------------------------------------------------------
!
        qp_temp = weight * cp * temp_curr
        call daxpy(cbs, qp_temp, BSCEval, 1, rhs, 1)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeLhsRigiTher(hhoCell, module_tang, weight, BSCEval, gbs, AT)
!
    implicit none
!
        type(HHO_Cell), intent(in)      :: hhoCell
        real(kind=8), intent(in)        :: module_tang(3,3)
        real(kind=8), intent(in)        :: weight
        real(kind=8), intent(in)        :: BSCEval(MSIZE_CELL_SCAL)
        integer, intent(in)             :: gbs
        real(kind=8), intent(inout)     :: AT(MSIZE_CELL_VEC, MSIZE_CELL_VEC)
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
        real(kind=8) :: qp_Agphi(MSIZE_CELL_VEC,3)
        integer:: i, k, gbs_cmp, col
! --------------------------------------------------------------------------------------------------
!
        gbs_cmp = gbs / hhoCell%ndim
! --------- Eval (A : gphi)_T
        call hhoComputeAgphi(hhoCell, module_tang, BSCEval, gbs, weight, qp_Agphi)
!
! -------- Compute scalar_product of (A_gphi, gphi)_T
        col = 1
        do i = 1, hhoCell%ndim
            do k = 1, gbs_cmp
                call daxpy(gbs, BSCEval(k), qp_Agphi(:,i), 1, AT(:, col), 1)
                col = col + 1
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
        real(kind=8), intent(in)        :: cp
        real(kind=8), intent(in)        :: weight
        real(kind=8), intent(in)        :: BSCEval(MSIZE_CELL_SCAL)
        integer, intent(in)             :: cbs
        real(kind=8), intent(inout)     :: lhs(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL)
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
! --------------------------------------------------------------------------------------------------
!
        coeff = cp * weight
        call dsyr('U', cbs, coeff, BSCEval, 1, lhs, MSIZE_TDOFS_SCAL)
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
        type(HHO_Cell), intent(in)      :: hhoCell
        integer, intent(in)             :: gbs
        real(kind=8), intent(in)        :: module_tang(3,3)
        real(kind=8), intent(in)        :: BSCEval(MSIZE_CELL_SCAL)
        real(kind=8), intent(in)        :: weight
        real(kind=8), intent(out)       :: Agphi(MSIZE_CELL_VEC,3)
!
! -----------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the scalar product (module_tang, gphi)_T

!   In hhoCell      : the current HHO Cell
!   In module_tang  : elasto_plastic moduli

!   In BSCEval      : Basis of one composant gphi
!   In gbs          : number of cols of Aphi
!   In weight       : quadrature weight
!   Out Agphi       : matriw of scalar product
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: qp_module_tang(3,3), qp_mod_vec(3)
        integer :: i, row, gbs_cmp, dim
! --------------------------------------------------------------------------------------------------
!
        Agphi = 0.d0
        dim = hhoCell%ndim
        gbs_cmp = gbs / dim
        qp_module_tang = weight * module_tang
!
        row = 1
        do i = 1, dim
            qp_mod_vec = qp_module_tang(:, i)
            call dger(gbs_cmp, dim, 1.d0, BSCEval, 1, qp_mod_vec, 1,&
                      Agphi(row:(row + gbs_cmp - 1), 1:dim), gbs_cmp)
            row = row + gbs_cmp
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function LambdaMax(fami, imate, npg, ndim, time) result(coeff)
!
    implicit none
!
        character(len=*), intent(in)  :: fami
        integer, intent(in)           :: imate, npg, ndim
        real(kind=8), intent(in)      :: time
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
                call rcvalb(fami, ipg, spt, '+', zi(imate),&
                            ' ', phenom, 1, 'INST', [time],&
                            1, nomres, valres, icodre, 1)
                lambda = valres(1)
                coeff = coeff + lambda
            else if (phenom.eq.'THER_ORTH') then
                ASSERT(ASTER_FALSE)
                nomres(1) = 'LAMBDA_L'
                nomres(2) = 'LAMBDA_T'
                nomres(3) = 'LAMBDA_N'
                call rcvalb(fami, ipg, spt, '+', zi(imate),&
                            ' ', phenom, 1, 'INST', [time],&
                            ndim, nomres, valres, icodre, 1)
                lambor(1) = valres(1)
                lambor(2) = valres(2)
                if(ndim == 3) lambor(3) = valres(3)
            else
                call utmess('F', 'ELEMENTS2_63')
            endif

        end do
!
        coeff = coeff / real(npg, kind=8)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCalcStabCoeff(hhoData, fami, nbQuadPoints, ndim, time)
!
    implicit none
!
        type(HHO_Data), intent(inout) :: hhoData
        character(len=4) :: fami
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
        if(hhoData%adapt()) then
            call jevech('PMATERC', 'L', imate)
            call hhoData%setCoeffStab( 10.d0*LambdaMax(fami, imate, nbQuadPoints, ndim, time))
       end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeBehaviourTher(phenom, fami, poum, kpg, ndim, time, imate, g_curr, &
                                       sig_curr, module_tang)
!
    implicit none
!
        character(len=32), intent(in) :: phenom
        character(len=8), intent(in) :: fami, poum
        integer, intent(in) :: imate, ndim, kpg
        real(kind=8), intent(in) :: time, g_curr(3)
        real(kind=8), intent(out) :: sig_curr(3), module_tang(3,3)
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - Integrate behaviour
!
! In hhoData          : information about the HHO formulation
! --------------------------------------------------------------------------------------------------
!
! --- Local variables
!
        character(len=16) :: nomres(3)
        real(kind=8) :: valres(3)
        integer :: icodre(3), idim
        integer, parameter :: spt = 1
        aster_logical :: aniso
        real(kind=8) :: lambda, lambor(3)
!
        sig_curr = 0.d0
        module_tang = 0.d0
!
! --- Eval lambda
!
        if (phenom .eq. 'THER') then
            nomres(1) = 'LAMBDA'
            call rcvalb(fami, kpg, spt, poum, zi(imate),&
                        ' ', phenom, 1, 'INST', [time],&
                        1, nomres, valres, icodre, 1)
            lambda = valres(1)
            aniso = ASTER_FALSE
        else if (phenom.eq.'THER_ORTH') then
            nomres(1) = 'LAMBDA_L'
            nomres(2) = 'LAMBDA_T'
            nomres(3) = 'LAMBDA_N'
            call rcvalb(fami, kpg, spt, poum, zi(imate),&
                        ' ', phenom, 1, 'INST', [time],&
                        ndim, nomres, valres, icodre, 1)
            lambor(1) = valres(1)
            lambor(2) = valres(2)
            if(ndim == 3) lambor(3) = valres(3)
            aniso = ASTER_TRUE
        else
            call utmess('F', 'ELEMENTS2_63')
        endif
!
! --- Rotation of local axis
!
        if(aniso) then
            ASSERT(ASTER_FALSE)
        else
            do idim = 1, ndim
                module_tang(idim, idim) = lambda
            end do
!
            sig_curr = - lambda * g_curr
        end if
    end subroutine
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoComputeRhoCpTher(phenom, fami, poum, kpg, time, imate, cp)
!
    implicit none
!
        character(len=32), intent(in) :: phenom
        character(len=8), intent(in) :: fami, poum
        integer, intent(in) :: imate,  kpg
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
            call rcvalb(fami, kpg, spt, poum, zi(imate),&
                        ' ', phenom, 1, 'INST', [time],&
                        1, 'RHO_CP', cp_, icod(1), 1)
        else if (phenom .eq. 'THER_ORTH') then
            call rcvalb(fami, kpg, spt, poum, zi(imate),&
                        ' ', phenom, 1, 'INST', [time],&
                        1, 'RHO_CP', cp_, icod(1), 1)
        else
            call utmess('F', 'ELEMENTS2_63')
        endif
!
        cp = cp_(1)
    end subroutine
!
end module
