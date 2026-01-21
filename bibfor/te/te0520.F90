! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine te0520(option, nomte)
!
    use coupling_penalisation_module
    use coupling_quadrature_module
    use coupling_type
    use FE_quadrature_module
    use FE_topo_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_penalisation_module.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/writeVector.h"
#include "FE_basis_module.h"
!
    character(len=16), intent(in) :: option, nomte
!
!===================================================================================================
!
!
! Mechanics - Coupling FEM/HHO with PENALISATION
!
!===================================================================================================
!
    type(FE_Skin):: feFace
    type(HHO_Face):: hhoFace
    type(HHO_Data):: hhoData
    type(FE_Quadrature) :: FEQuad
    type(HHO_Quadrature) :: hhoQuad
    type(CouplingMap) :: cplPenaMap
    type(CouplingData) :: cplPenaData
    real(kind=8) :: rhs(MSIZE_CPL_PENA)
    type(HHO_matrix) :: lhs
!
! --- Initialize topologie
!
    call cplPenaInitTopo(feFace, hhoFace, hhoData)
!
! --- Initialize mapping
!
    call cplPenaInitMap(feFace, hhoFace, hhoData, cplPenaMap)
!
! --- Initialize quadrature
!
    call cplGetQuadFEMHHO(feFace, hhoFace, hhoData, FEQuad, hhoQuad)
!
! --- Initialize data
!
    call cplPenaInitData(option, cplPenaMap, cplPenaData)
!
! - Computation
!
    if (option == "CHAR_MECA_CPL") then
!
! --- Compute coupling residual
!
        call cplPenaRhs(feFace, hhoFace, hhoData, FEQuad, hhoQuad, &
                        cplPenaData, rhs)
!
! --- Write vector
!
        call writeVector('PVECTUR', cplPenaMap%nbDoFs, rhs)
!
    elseif (option == "RIGI_CPL") then
!
! --- Compute coupling matrix
!
        call cplPenaLhs(feFace, hhoFace, hhoData, FEQuad, hhoQuad, &
                        cplPenaData, lhs)
!
! - Write matrix
!
        call lhs%write("PMATUUR", ASTER_TRUE)
        call lhs%free()
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
