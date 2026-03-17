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
subroutine te0495(nomopt, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use HHO_type
    use HHO_compor_module
    use HHO_utils_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_Meca_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_LargeStrainMeca_module
    use HHO_basis_module
    use HHO_eval_module
    use HHO_matrix_module
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/pidefo.h"
#include "asterfort/pielas.h"
#include "asterfort/readVector.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - STAT_NON_LINE - Pilotage
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
    type(HHO_Quadrature) :: hhoQuadCellRigi
    type(HHO_Compor_State) :: hhoCS
    type(HHO_basis_cell) :: hhoBasisCell
    type(Behaviour_Integ) :: BEHinteg
!
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
    real(kind=8), dimension(MSIZE_TDOFS_VEC) :: depl_0r, depl_1r
    real(kind=8), dimension(MSIZE_CELL_MAT) :: E_prev_coeff, E_incr_coeff
    real(kind=8), dimension(MSIZE_CELL_MAT) :: E_0r_coeff, E_1r_coeff
    real(kind=8) :: tau, etamin, etamax, time_prev, time_pilo
    real(kind=8) :: coorpg(3), E_prev(6), E_pilo(6), E_1(6), E_0(6), E_incr(6)
    real(kind=8) :: copilo(5, MAX_QP_CELL), sigma(6)
    real(kind=8) :: G_prev(3, 3), G_incr(3, 3), G_1(3, 3), G_0(3, 3)
    real(kind=8) :: F_prev(3, 3), F_incr(3, 3), F_1(3, 3), F_0(3, 3)
    integer(kind=8) :: cbs, fbs, total_dofs, gbs, gbs_sym
    integer(kind=8) :: ipg, npg, k
    integer(kind=8) :: iborne, ictau, itype, imate
    character(len=4), parameter :: fami = 'RIGI'
    character(len=16) :: pilo
    real(kind=8), pointer :: v_copilo(:) => null()
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
! --- Get element parameters
!
    call elrefe_info(fami=fami, npg=npg)
!
! --- Number of dofs
    call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs, gbs_sym)
    ASSERT(total_dofs <= MSIZE_TDOFS_VEC)
!
    if (nomopt /= "PILO_PRED_DEFO" .and. nomopt /= "PILO_PRED_ELAS") then
        ASSERT(ASTER_FALSE)
    end if
!
! --- Initialize quadrature for the rigidity
!
    call hhoQuadCellRigi%initCell(hhoCell, npg)
!
! --- Type of finite element
!
    call hhoCS%initialize(fami, nomopt, hhoCell%ndim, hhoCell%barycenter)
    call hhoMecaState%initialize(hhoCell, hhoData, hhoCS)
    call hhoBasisCell%initialize(hhoCell)
!
! --- Compute Operators
!
    if (hhoData%precompute()) then
        call hhoReloadPreCalcMeca(hhoCell, hhoData, hhoCS%l_largestrain, hhoMecaState%grad)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, hhoCS%l_largestrain, hhoMecaState%grad)
    end if
!
    call readVector('PDEPL0R', total_dofs, depl_0r)
    call readVector('PDEPL1R', total_dofs, depl_1r)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PTYPEPI', 'L', itype)
!
    pilo = zk16(itype)
    if (pilo .eq. 'PRED_ELAS') then
        call jevech('PCDTAU', 'L', ictau)
        call jevech('PBORNPI', 'L', iborne)
!
        tau = zr(ictau)
        etamin = zr(iborne+1)
        etamax = zr(iborne)
    end if
!
! - Continuation method: no time !
    time_prev = r8vide()
    time_pilo = r8vide()
!
! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)
!
! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(hhoCell%ndim, hhoCS%typmod, hhoCS%option, &
                              hhoCS%compor, hhoCS%carcri, &
                              time_prev, time_pilo, &
                              hhoCS%fami, hhoCS%imater, &
                              BEHinteg)
!
! - Prepare external state variables (geometry)
    call behaviourPrepESVAGeomHHO(hhoCell, hhoQuadCellRigi, BEHinteg)
!
! ----- compute E_prev = gradrec_sym * depl_prev
    call hhoMecaState%grad%dot(hhoMecaState%depl_prev, E_prev_coeff)
    call hhoMecaState%grad%dot(hhoMecaState%depl_incr, E_incr_coeff)
    call hhoMecaState%grad%dot(depl_0r, E_0r_coeff)
    call hhoMecaState%grad%dot(depl_1r, E_1r_coeff)
!
    call matfpe(-1)
    copilo = r8vide()
    sigma = 0.d0
!
    do ipg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
! --------- Eval basis function at the quadrature point
        call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
        if (hhoCS%l_largestrain) then
            G_prev = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, E_prev_coeff)
            call hhoCalculF(G_prev, F_prev)
            call hhoCalculGreenLagrange(hhoCell%ndim, F_prev, E_prev)
!
            G_incr = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, E_incr_coeff)
            call hhoCalculF(G_incr, F_incr)
            call hhoCalculGreenLagrange(hhoCell%ndim, F_incr, E_incr)
!
            G_0 = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, E_0r_coeff)
            call hhoCalculF(G_0, F_0)
            call hhoCalculGreenLagrange(hhoCell%ndim, F_0, E_0)
!
            E_pilo = E_incr+E_0
!
            G_1 = hhoEvalMatCell(hhoCell%ndim, gbs, BSCEval, E_1r_coeff)
            call hhoCalculF(G_1, F_1)
            call hhoCalculGreenLagrange(hhoCell%ndim, F_1, E_1)
        else
! --------- Eval deformations
            E_prev = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_prev_coeff)
            E_incr = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_incr_coeff)
            E_0 = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_0r_coeff)
            E_pilo = E_incr+E_0
            E_1 = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_1r_coeff)
        end if
!
! --- PILOTAGE PAR L'INCREMENT DE DEFORMATION
!
        if (pilo .eq. 'DEFORMATION') then
!
            call pidefo(hhoCell%ndim, npg, ipg, &
                        hhoCS%compor, F_prev, &
                        E_prev, E_pilo, E_1, copilo)
!
! --- PILOTAGE PAR LA PREDICTION ELASTIQUE
!
        else if (pilo .eq. 'PRED_ELAS') then
            sigma(1:hhoCS%nbsigm) = hhoCS%sig_prev((ipg-1)*hhoCS%nbsigm+1:ipg*hhoCS%nbsigm)
            do k = 4, hhoCS%nbsigm
                sigma(k) = sigma(k)*rac2
            end do
            call pielas(BEHinteg, hhoCell%ndim, npg, ipg, hhoCS%compor, &
                        hhoCS%typmod, zi(imate), hhoCS%lgpg, hhoCS%vari_prev, &
                        E_prev, E_pilo, E_1, sigma, etamin, etamax, &
                        tau, copilo)
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
    call matfpe(1)
!
    call jevech('PCOPILO', 'E', vr=v_copilo)
!
    do ipg = 1, npg
        v_copilo((ipg-1)*5+1:ipg*5) = copilo(1:5, ipg)
    end do
!
end subroutine
