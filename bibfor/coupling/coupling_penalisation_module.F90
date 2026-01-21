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
module coupling_penalisation_module
!
    use coupling_computation_module
    use coupling_quadrature_module
    use coupling_type
    use FE_algebra_module
    use FE_Basis_module
    use FE_quadrature_module
    use FE_topo_module
    use HHO_init_module
    use HHO_matrix_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
!
    implicit none
!
    private
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/coupling_penalisation_module.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/readVector.h"
#include "FE_basis_module.h"
!
! --------------------------------------------------------------------------------------------------
!
! Coupling - PENALISATION methods
!
! --------------------------------------------------------------------------------------------------
!
!
!
    public :: cplPenaInitTopo, cplPenaInitMap, cplPenaInitData
    public :: cplPenaRhs, cplPenaLhs
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaInitTopo(feFace, hhoFace, hhoData)
!
        implicit none
!
        type(FE_Skin), intent(out) :: feFace
        type(HHO_Face), intent(out) :: hhoFace
        type(HHO_Data), intent(out) :: hhoData
!
!===================================================================================================
!    Initialize FEM and HHO topology - Penalisation
!===================================================================================================
!
        call feFace%init()
!
        call hhoInfoInitFace(hhoFace, hhoData)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaInitMap(FEFace, hhoFace, hhoData, cplMap)
!
        implicit none
!
        type(CouplingMap), intent(inout) :: cplMap
        type(FE_Skin), intent(in) :: FEFace
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Data), intent(in) :: hhoData
!
!===================================================================================================
!    Initialize FEM and HHO topology
!===================================================================================================
!
        type(FE_basis) :: FEBasis
        integer(kind=8) :: iBase, idof, idim, fbs
        call FEBasis%initFace(FEFace)
        ASSERT(FEBasis%typeEF == EF_LAGRANGE)
!
        idof = 1
!
        ASSERT(FEBasis%size == FEFace%nbnodes)
        do iBase = 1, FEBasis%size
            do idim = 1, FEBasis%ndim
                cplMap%mapDoFsFEFace(idof) = idof
                idof = idof+1
            end do
        end do
!
        call hhoMecaFaceDofs(hhoFace, hhoData, fbs)
!
        do iBase = 1, fbs
            cplMap%mapDoFshhoFace(iBase) = idof
            idof = idof+1
        end do
!
        cplMap%nbDoFs = FEBasis%size*FEBasis%ndim+fbs
        cplMap%nbDoFsFECell = 0
        cplMap%nbDoFsFEFace = FEBasis%size*FEBasis%ndim
        cplMap%nbDoFshhoFace = fbs
        ASSERT(cplMap%nbDoFs <= MSIZE_CPL_PENA)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaInitData(option, cplMap, cplData)
!
        implicit none
!
        character(len=16), intent(in) :: option
        type(CouplingMap), intent(in) :: cplMap
        type(CouplingData), intent(inout) :: cplData
!
!===================================================================================================
!    Initialize FEM and HHO data
!===================================================================================================
!
        integer(kind=8) :: jpair
!
        if (option == "CHAR_MECA_CPL") then
            call readVector('PDEPLMR', cplMap%nbDoFs, cplData%disp_prev)
            call readVector('PDEPLPR', cplMap%nbDoFs, cplData%disp_curr)
            call daxpy_1(cplMap%nbDoFs, 1.d0, cplData%disp_prev, cplData%disp_curr)
        end if
!
        call jevech('PPAIRR', 'L', jpair)
        cplData%coef_pena = zr(jpair-1+30)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaRhs(feFace, hhoFace, hhoData, FEQuad, hhoQuad, cplData, rhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: feFace
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Data), intent(in) :: hhoData
        type(FE_Quadrature), intent(in) :: FEQuad
        type(HHO_Quadrature), intent(in) :: hhoQuad
        type(CouplingData), intent(in) :: cplData
        real(kind=8), intent(inout) :: rhs(MSIZE_CPL_PENA)
!
!===================================================================================================
!    PENALISATION - FEM/HHO - Compute residual
!===================================================================================================
!
        type(HHO_matrix) :: mat
!
        rhs = 0.d0
!
        call cplCmpPenaFEHHOMatVec(feFace, hhoFace, hhoData, FEQuad, hhoQuad, &
                                   cplData%coef_pena, mat)
        call mat%dot(cplData%disp_curr, rhs)
        call mat%free()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
!
    subroutine cplPenaLhs(feFace, hhoFace, hhoData, FEQuad, hhoQuad, cplData, lhs)
!
        implicit none
!
        type(FE_Skin), intent(in) :: feFace
        type(HHO_Face), intent(in) :: hhoFace
        type(HHO_Data), intent(in) :: hhoData
        type(FE_Quadrature), intent(in) :: FEQuad
        type(HHO_Quadrature), intent(in) :: hhoQuad
        type(CouplingData), intent(in) :: cplData
        type(HHO_matrix), intent(inout) :: lhs
!
!===================================================================================================
!    PENALISATION - FEM/HHO - Compute matrix
!===================================================================================================
!
        call cplCmpPenaFEHHOMatVec(feFace, hhoFace, hhoData, FEQuad, hhoQuad, &
                                   cplData%coef_pena, lhs)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
end module
