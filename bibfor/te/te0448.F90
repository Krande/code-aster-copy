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
subroutine te0448(nomopt, nomte)
!
    use HHO_basis_module
    use HHO_eval_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_gradrec_module
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use HHO_utils_module
    use HHO_matrix_module
    use HHO_algebra_module
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
#include "asterfort/readVector.h"
#include "blas/dgemv.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - EPSI_ELGA
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
    integer(kind=8) :: cbs, fbs, total_dofs, gbs, gbs_sym
    integer(kind=8) :: npg
    integer(kind=8) :: ipg, idefo, nsig
    aster_logical :: l_largestrains
    character(len=4) :: fami
    character(len=8) :: typmod(2)
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    real(kind=8) :: G_curr(3, 3), E_curr(6)
    real(kind=8) :: coorpg(3)
    real(kind=8) :: BSCEval(MSIZE_CELL_SCAL)
    real(kind=8), dimension(MSIZE_TDOFS_VEC) :: depl_curr
    real(kind=8), dimension(MSIZE_CELL_MAT) :: G_curr_coeff
    type(HHO_matrix) :: gradrec
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
! --- Get element parameters
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, npg=npg)
!
! --- Number of dofs
    call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, &
                       gbs, gbs_sym)
    nsig = nbsigm()
    ASSERT(cbs <= MSIZE_CELL_VEC)
    ASSERT(fbs <= MSIZE_FACE_VEC)
    ASSERT(total_dofs <= MSIZE_TDOFS_VEC)
!
    ASSERT(nomopt .eq. 'EPSI_ELGA')
!
! --- Initialize quadrature for the rigidity
    call hhoQuadCellRigi%initCell(hhoCell, npg)
!
! --- Type of finite element
!
    select case (hhoCell%ndim)
    case (3)
        typmod(1) = '3D'
    case (2)
        if (lteatt('AXIS', 'OUI')) then
            ASSERT(ASTER_FALSE)
            typmod(1) = 'AXIS'
        else if (lteatt('C_PLAN', 'OUI')) then
            ASSERT(ASTER_FALSE)
            typmod(1) = 'C_PLAN'
        else if (lteatt('D_PLAN', 'OUI')) then
            typmod(1) = 'D_PLAN'
        else
            ASSERT(ASTER_FALSE)
        end if
    case default
        ASSERT(ASTER_FALSE)
    end select
    typmod(2) = 'HHO'
!
    call jevech('PDEFOPG', 'E', idefo)
!
! --- Large strains ?
!
    l_largestrains = ASTER_FALSE
!
! --- Compute Operators
!
    if (l_largestrains) then
!
! ----- Compute Gradient reconstruction
        call hhoGradRecFullMat(hhoCell, hhoData, gradrec)
    else
!
! ----- Compute Symmetric Gradient reconstruction
        call hhoGradRecSymFullMat(hhoCell, hhoData, gradrec)
    end if
!
! --- get displacement
!
    depl_curr = 0.d0
    call readVector('PDEPLAR', total_dofs, depl_curr)
!
! ----- init basis
!
    call hhoBasisCell%initialize(hhoCell)
!
! --- Compute local contribution
!
    call hho_dgemv_N(1.d0, gradrec, depl_curr, 0.d0, G_curr_coeff)
!
! ----- Loop on quadrature point
!
    do ipg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, ipg)
!
! --------- Eval basis function at the quadrature point
!
        call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
        if (l_largestrains) then
            G_curr = hhoEvalMatCell( &
                     hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), G_curr_coeff, gbs)
        else
            E_curr = hhoEvalSymMatCell( &
                     hhoBasisCell, hhoData%grad_degree(), coorpg(1:3), G_curr_coeff, gbs_sym)
            zr(idefo-1+(ipg-1)*nsig+1:idefo-1+ipg*nsig) = E_curr(1:nsig)
        end if
    end do
!
    call gradrec%free()
!
end subroutine
