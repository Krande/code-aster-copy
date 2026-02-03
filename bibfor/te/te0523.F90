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
subroutine te0523(option, nomte)
!
    use coupling_lagrangian_module
    use coupling_quadrature_module
    use coupling_type
    use FE_quadrature_module
    use FE_topo_module
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_penalisation_module.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/writeVector.h"
#include "FE_basis_module.h"
!
    character(len=16), intent(in) :: option, nomte
!
!===================================================================================================
!
!
! Mechanics - Coupling FEM/FEM with augmented Lagrangian
!
!===================================================================================================
!
    type(FE_Skin):: FEFaceSl, FEFaceMa, FEFaceLagSl
    type(FE_Quadrature) :: FEQuadSl, FEQuadMa

    type(CouplingMap) :: cplLagrMap
    type(CouplingData) :: cplLagrData
    real(kind=8) :: rhs(MSIZE_CPL_LAGR_FEM)
    real(kind=8) :: lhs(MSIZE_CPL_LAGR_FEM, MSIZE_CPL_LAGR_FEM)
!
! --- Initialize topologie
!
    call cplLagrInitTopoFEMFEM(FEFaceSl, FEFaceMa, FEFaceLagSl)
!
! --- Initialize mapping
!
    call cplLagrInitMapFEMFEM(FEFaceSl, FEFaceMa, FEFaceLagSl, cplLagrMap)
!
! --- Initialize quadrature
!
    call cplGetQuadFEMFEM(FEFaceSl, FEFaceMa, FEQuadSl, FEQuadMa)
!
! --- Initialize data
!
    call cplLagrInitData(option, cplLagrMap, cplLagrData)
!
! - Computation
!
    if (option == "CHAR_MECA_CPL") then
!
! --- Compute coupling residual
!
        call cplLagrRhsFEMFEM(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, &
                              cplLagrMap, cplLagrData, rhs)
!
! --- Write vector
!
        call writeVector('PVECTUR', cplLagrMap%nbDoFs, rhs)
!
    elseif (option == "RIGI_CPL" .or. option == "RIGI_ELAS_CPL") then
!
! --- Compute coupling matrix
!
        call cplLagrLhsFEMFEM(FEFaceSl, FEFaceMa, FEFaceLagSl, FEQuadSl, FEQuadMa, &
                              cplLagrMap, cplLagrData, lhs)
!
! - Write matrix
!
        call writeMatrix("PMATUUR", cplLagrMap%nbDoFs, cplLagrMap%nbDoFs, &
                         ASTER_TRUE, lhs)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
