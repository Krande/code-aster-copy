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

subroutine te0450(nomopt, nomte)
!
    use HHO_basis_module
    use HHO_compor_module
    use HHO_eval_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
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
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - FORC_NODA and REFE_FORC_NODA
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Quadrature) :: hhoQuadCellRigi
    type(HHO_Compor_State) :: hhoCS
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
!
    integer :: cbs, fbs, total_dofs, npg
    aster_logical :: l_largestrains
    character(len=4), parameter :: fami = "RIGI"
    real(kind=8) :: rhs(MSIZE_TDOFS_VEC)
!
! --- Get HHO informations
!
    call elrefe_info(fami=fami, npg=npg)
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCellRigi)
!
! --- Get element parameters
!
    rhs = 0.d0
!
! --- Number of dofs
    call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
! --- Type of finite element
!
    call hhoCS%initialize(fami, nomopt, hhoCell%ndim, hhoCell%barycenter)
    call hhoMecaState%initialize(hhoCell, hhoData, hhoCS)
!
! --- Large strains ?
!
    l_largestrains = hhoCS%l_largestrain
!
! --- Compute Operators
!
    if (hhoData%precompute()) then
        call hhoReloadPreCalcMeca(hhoCell, hhoData, l_largestrains, &
                                  hhoMecaState%grad, hhoMecaState%stab)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, l_largestrains, hhoMecaState%grad, hhoMecaState%stab)
    end if
!
    if (nomopt == "FORC_NODA") then
        call hhoLocalForcNoda(hhoCell, hhoData, hhoQuadCellRigi, hhoMecaState, &
                              hhoCS, hhoCS%sig_prev, rhs)
    elseif (nomopt == "REFE_FORC_NODA") then
        ASSERT(ASTER_FALSE)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! --- Save rhs
!
    call hhoRenumMecaVec(hhoCell, hhoData, rhs)
    call writeVector('PVECTUR', total_dofs, rhs)
!
end subroutine
