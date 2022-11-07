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
subroutine te0473(option, nomte)
!
use HHO_type
use HHO_size_module
use HHO_init_module, only : hhoInfoInitCell
use HHO_L2proj_module
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Thermics - HHO_PROJ_THER
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16), intent(in) :: option, nomte
!
! --- Local variables
!
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    integer :: cbs, fbs, total_dofs, jvale
    real(kind=8) :: value
    real(kind=8), dimension(MSIZE_TDOFS_SCAL) :: coeff_L2Proj
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
    call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
    call jevech("PVALE_R", "L", jvale)
    value = zr(jvale)
!
    call hhoL2ProjScal(hhoCell, hhoData, value, coeff_L2Proj)
!
    call writeVector("PTEMP_R", total_dofs, coeff_L2Proj)
!
end subroutine
