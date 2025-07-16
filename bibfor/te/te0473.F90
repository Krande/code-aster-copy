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
subroutine te0473(option, nomte)
!
    use HHO_type
    use HHO_size_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_L2proj_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/writeVector.h"
#include "asterfort/assert.h"
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
    integer(kind=8) :: cbs, fbs, total_dofs, jvale, jfunc
    real(kind=8) :: time
    character(len=8) :: func
    real(kind=8), dimension(MSIZE_TDOFS_VEC) :: coeff_L2Proj
    aster_logical :: with_faces
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
    call jevech("PFUNC_R", "L", jfunc)
    call jevech("PINSTPR", "L", jvale)
    time = zr(jvale)
!
    if (option == "HHO_PROJ_THER") then
!
        func = zk8(jfunc-1+1)
        with_faces = zk8(jfunc-1+2) == "ALL"
        call hhoL2ProjScal(hhoCell, hhoData, func, time, coeff_L2Proj, with_faces)
!
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        call writeVector("PTEMP_R", total_dofs, coeff_L2Proj)
    elseif (option == "HHO_PROJ_MECA") then
        with_faces = zk8(jfunc-1+hhoCell%ndim+1) == "ALL"
        call hhoL2ProjVec(hhoCell, hhoData, zk8(jfunc), time, coeff_L2Proj, with_faces)
!
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        call writeVector("PDEPL_R", total_dofs, coeff_L2Proj)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
