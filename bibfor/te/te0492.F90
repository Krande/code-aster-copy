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

subroutine te0492(nomopt, nomte)
!
    use HHO_type
    use HHO_size_module, only: hhoTherDofs
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_Dirichlet_module
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/writeVector.h"
#include "asterfort/HHO_size_module.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO - Thermics
!  Option: AFFE_CHAR_CINE_R
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: nomte, nomopt
!
! -- Local variables

    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    real(kind=8) :: rhs_cine(MSIZE_TDOFS_SCAL)
    integer(kind=8) :: cbs, fbs, total_dofs
    real(kind=8), pointer :: r_vale(:) => null()
!
! --- Retrieve HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
    if (nomopt .eq. 'HHO_CINE_R_THER') then
!
! --- Read Name of field
!
        call jevech('PCMPVALE', 'L', vr=r_vale)
!
! --- Projection of the boundary conditions
!
        call hhoDiriTherProjReal(hhoCell, hhoData, r_vale, rhs_cine)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
! -- Save
!
    call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    call writeVector('PCINE', total_dofs, rhs_cine)
!
end subroutine
