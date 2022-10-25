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

subroutine te0454(nomopt, nomte)
!
use HHO_type
use HHO_utils_module
use HHO_size_module
use HHO_quadrature_module
use HHO_ther_module
use HHO_init_module, only : hhoInfoInitCell
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/writeMatrix.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - RIGI_THER
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
    character(len=4) :: fami
    character(len=8) :: typmod(2)
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    real(kind=8), dimension(MSIZE_CELL_VEC, MSIZE_TDOFS_SCAL)   :: gradfull
    real(kind=8), dimension(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL)  :: lhs, stab
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
    call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    ASSERT(cbs <= MSIZE_CELL_SCAL)
    ASSERT(fbs <= MSIZE_FACE_SCAL)
    ASSERT(total_dofs <= MSIZE_TDOFS_SCAL)
!
    if (nomopt /= "RIGI_THER") then
        ASSERT(ASTER_FALSE)
    end if
!
! --- Initialize quadrature for the rigidity
!
    call hhoQuadCellRigi%initCell(hhoCell, npg)
!
! --- Type of finite element
!
    select case (hhoCell%ndim)
        case(3)
            typmod(1) = '3D'
        case (2)
            typmod(1) = 'PLAN'
        case default
            ASSERT(ASTER_FALSE)
    end select
    typmod(2) = 'HHO'
!
! --- Compute Operators
!
    call hhoCalcOpTher(hhoCell, hhoData, gradfull, stab)
!
! --- Compute local contribution
!
    call hhoLocalRigiTher(hhoCell, hhoData, hhoQuadCellRigi, gradfull, stab, &
                             fami, lhs=lhs)
!
! --- Save lhs
!
    call hhoRenumTherMat(hhoCell, hhoData, lhs)
    call writeMatrix('PMATTTR', total_dofs, total_dofs, ASTER_TRUE, lhs)
!
end subroutine
