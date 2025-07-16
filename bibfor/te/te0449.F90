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
subroutine te0449(nomopt, nomte)
!
    use HHO_type
    use HHO_utils_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_ther_module
    use HHO_basis_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_matrix_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/lteatt.h"
#include "jeveux.h"
!
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - MASS_THER
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Quadrature) :: hhoQuadCellMass
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    integer(kind=8) :: npg
    character(len=8), parameter :: fami = 'MASS'
    type(HHO_matrix) :: lhs
!
! --- Get element parameters
!
    call elrefe_info(fami=fami, npg=npg)
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCellMass)
!
    if (nomopt /= "MASS_THER") then
        ASSERT(ASTER_FALSE)
    end if
!
    if (lteatt('LUMPE', 'OUI')) then
        ASSERT(ASTER_FALSE)
    end if
!
    call hhoLocalMassTher(hhoCell, hhoData, hhoQuadCellMass, fami, lhs=lhs)
!
! --- Save lhs
!
    call lhs%write('PMATTTR', ASTER_TRUE)
    call lhs%free()
!
end subroutine
