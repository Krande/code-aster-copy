! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine te0477(nomopt, nomte)
!
    use HHO_type
    use HHO_utils_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_ther_module
    use HHO_init_module, only: hhoInfoInitCell
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/writeVector.h"
#include "asterfort/teattr.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - RAPH_THER
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Quadrature) :: hhoQuadCellRigi
    integer :: cbs, fbs, total_dofs, npg
    character(len=8), parameter :: fami = 'RIGI'
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    real(kind=8), dimension(MSIZE_CELL_VEC, MSIZE_TDOFS_SCAL) :: gradfull
    real(kind=8), dimension(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL) :: stab
    real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: rhs
!
! --- Get element parameters
!
    call elrefe_info(fami=fami, npg=npg)
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCellRigi)
!
! --- Number of dofs
    call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    ASSERT(cbs <= MSIZE_CELL_SCAL)
    ASSERT(fbs <= MSIZE_FACE_SCAL)
    ASSERT(total_dofs <= MSIZE_TDOFS_SCAL)
!
    if (nomopt /= "RAPH_THER") then
        ASSERT(ASTER_FALSE)
    end if
!
! --- Compute Operators
!
    call hhoCalcOpTher(hhoCell, hhoData, gradfull, stab)
!
! --- Compute local contribution
!
    call hhoLocalRigiTher(hhoCell, hhoData, hhoQuadCellRigi, nomopt, gradfull, stab, &
                          fami, rhs=rhs)
!
! --- Save lhs
!
    call hhoRenumTherVec(hhoCell, hhoData, rhs)
    call writeVector('PRESIDU', total_dofs, rhs)
!
end subroutine
