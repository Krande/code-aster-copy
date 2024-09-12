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
!
subroutine te0450(nomopt, nomte)
!
    use HHO_basis_module
    use HHO_compor_module
    use HHO_eval_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_LargeStrainMeca_module, only: hhoComputeRhsLarge, hhoCalculF
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_SmallStrainMeca_module, only: hhoComputeRhsSmall, tranfoMatToSym
    use HHO_type
    use HHO_utils_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sigtopk1.h"
#include "asterfort/writeVector.h"
#include "asterfort/readVector.h"
#include "blas/dcopy.h"
#include "blas/dgemm.h"
#include "blas/dgemv.h"
#include "blas/dscal.h"
#include "blas/dsymv.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - FORC_NODA
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_basis_cell) :: hhoBasisCell
    type(HHO_Quadrature) :: hhoQuadCellRigi
    type(HHO_Compor_State) :: hhoCS
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
!
    integer :: cbs, fbs, total_dofs, gbs, gbs_sym
    integer :: npg, gbs_curr, gbs_cmp
    integer :: ipg, ncomp
    aster_logical :: l_largestrains
    character(len=4) :: fami
    real(kind=8) :: Cauchy_curr(6), PK1_curr(3, 3), G_curr(3, 3), F_curr(3, 3)
    real(kind=8) :: rhs(MSIZE_TDOFS_VEC), coorpg(3), weight
    real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
    real(kind=8), dimension(MSIZE_CELL_MAT) :: bT, G_curr_coeff
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
! --- Get element parameters
!
    Cauchy_curr = 0.d0
    bT = 0.d0
    rhs = 0.d0
    fami = 'RIGI'
    call elrefe_info(fami=fami, npg=npg)
!
! --- Number of dofs
    call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                       gbs, gbs_sym)
    gbs_cmp = gbs/(hhoCell%ndim*hhoCell%ndim)
    ASSERT(cbs <= MSIZE_CELL_VEC)
    ASSERT(fbs <= MSIZE_FACE_VEC)
    ASSERT(total_dofs <= MSIZE_TDOFS_VEC)
!
    ASSERT(nomopt .eq. 'FORC_NODA')
!
! --- Initialize quadrature for the rigidity
    call hhoQuadCellRigi%initCell(hhoCell, npg)
!
! --- Type of finite element
!
    call hhoCS%initialize(fami, nomopt, hhoCell%ndim, hhoCell%barycenter)
    call hhoMecaState%initialize(hhoCell, hhoData, hhoCS)
!
! --- Properties of behaviour
!
    ncomp = hhoCS%nbsigm
!
! --- Large strains ?
!
    l_largestrains = hhoCS%l_largestrain
!
! --- Compute Operators
!
    call hhoCalcOpMeca(hhoCell, hhoData, l_largestrains, hhoMecaState%grad, hhoMecaState%stab)
!
! ----- init basis
!
    call hhoBasisCell%initialize(hhoCell)
!
! --- Compute local contribution
!
    if (l_largestrains) then
        call dgemv('N', gbs, total_dofs, 1.d0, hhoMecaState%grad, &
                   MSIZE_CELL_MAT, hhoMecaState%depl_curr, 1, 0.d0, G_curr_coeff, &
                   1)
        gbs_curr = gbs
    else
        gbs_curr = gbs_sym
    end if
!
! ----- Loop on quadrature point
!
    do ipg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
        weight = hhoQuadCellRigi%weights(ipg)
!
! -------- tranform sigm in symmetric form
!
        call tranfoMatToSym(hhoCell%ndim, hhoCS%sig_prev((ipg-1)*ncomp+1:ipg*ncomp), Cauchy_curr)
!
! --------- Eval basis function at the quadrature point
!
        call hhoBasisCell%BSEval(hhoCell, coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
! --------- Eval gradient at T- and T+
!
!
        if (l_largestrains) then
            G_curr = hhoEvalMatCell( &
                     hhoCell, hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), G_curr_coeff, gbs &
                     )
!
! --------- Eval gradient of the deformation at T- and T+
!
            call hhoCalculF(hhoCell%ndim, G_curr, F_curr)
!
            call sigtopk1(hhoCell%ndim, Cauchy_curr, F_curr, PK1_curr)
!
            call hhoComputeRhsLarge(hhoCell, PK1_curr, weight, BSCEval, gbs, &
                                    bT)
        else
!
            call hhoComputeRhsSmall(hhoCell, Cauchy_curr, weight, BSCEval, gbs_cmp, &
                                    bT)
        end if
    end do
!
    call dgemv('T', gbs_curr, total_dofs, 1.d0, hhoMecaState%grad, &
               MSIZE_CELL_MAT, bT, 1, 1.d0, rhs, &
               1)
!
! --- add stabilization
!
    call hhoCalcStabCoeffMeca(hhoData, hhoCS%fami, 0.d0, hhoQuadCellRigi)
!
    call dsymv('U', total_dofs, hhoData%coeff_stab(), hhoMecaState%stab, MSIZE_TDOFS_VEC, &
               hhoMecaState%depl_curr, 1, 1.d0, rhs, 1)
!
! --- Save rhs
!
    call hhoRenumMecaVec(hhoCell, hhoData, rhs)
    call writeVector('PVECTUR', total_dofs, rhs)
!
end subroutine
