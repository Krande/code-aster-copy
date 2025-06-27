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
subroutine te0052(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_eval_module
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"

    character(len=16) :: option, nomte
!
!     BUT:
!       CALCUL DE GRADIENT DE TEMPERATURE AUX POINTS DE GAUSS
!       OPTION : 'GRAT_ELGA'
!
! ---------------------------------------------------------------------
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) ::kp
    real(kind=8) :: dtpg(3)
    real(kind=8), pointer :: tempi(:) => null()
    real(kind=8), pointer :: gradt(:) => null()
!
! ----------------------------------------------------------------------
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
!
    call jevech('PGRATPG', 'E', vr=gradt)
    call jevech('PTEMPER', 'L', vr=tempi)
!
    do kp = 1, FEQuadCell%nbQuadPoints
        dtpg = FEEvalGradVec(FEBasis, tempi, FEQuadCell%points_param(1:3, kp))
        gradt(FECell%ndim*(kp-1)+1:FECell%ndim*(kp-1)+FECell%ndim) = dtpg(1:FECell%ndim)
    end do
!
end subroutine
