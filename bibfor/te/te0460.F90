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

subroutine te0460(nomopt, nomte)
!
    use HHO_type
    use HHO_size_module
    use HHO_stabilization_module, only: hhoStabScal, hdgStabScal
    use HHO_gradrec_module, only: hhoGradRecVec, hhoGradRecFullVec
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_matrix_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/assert.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - Precomputation of operators
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    integer(kind=8) :: cbs, fbs, total_dofs, gbs
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_matrix) :: gradfullvec, stabscal, gradrec_scal

    ASSERT(nomopt == 'HHO_PRECALC_OP')
!
! --- Retrieve HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
    call hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)

! ----- Compute vectoriel Gradient reconstruction
    call hhoGradRecFullVec(hhoCell, hhoData, gradfullvec)
!
    call gradfullvec%write('PCHHOGT', ASTER_FALSE)
    call gradfullvec%free()
!
! ----- Compute Stabilizatiion
!
    if (hhoData%cell_degree() <= hhoData%face_degree()) then
        call hhoGradRecVec(hhoCell, hhoData, gradrec_scal)
        call hhoStabScal(hhoCell, hhoData, gradrec_scal, stabscal)
        call gradrec_scal%free()
    else if (hhoData%cell_degree() == (hhoData%face_degree()+1)) then
        call hdgStabScal(hhoCell, hhoData, stabscal)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call stabscal%write('PCHHOST', ASTER_TRUE)
    call stabscal%free()
!
end subroutine
