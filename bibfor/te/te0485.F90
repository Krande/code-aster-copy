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

subroutine te0485(nomopt, nomte)
!
    use Behaviour_module, only: behaviourOption
    use HHO_type
    use HHO_utils_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_Meca_module
    use HHO_compor_module
    use HHO_GV_module
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
#include "asterfort/jevech.h"
#include "asterfort/nmtstm.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!  HHO
!  Mechanics - STAT_NON_LINE - GRAD_VARI
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! --------------------------------------------------------------------------------------------------
    character(len=16) :: nomte, nomopt
!
! --- Local variables
!
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
    type(HHO_GV_State) :: hhoGVState
    type(HHO_Compor_State) :: hhoComporState
    type(HHO_Quadrature) :: hhoQuadCellRigi
    integer(kind=8) :: mk_cbs, mk_fbs, mk_total_dofs
    integer(kind=8) :: gv_cbs, gv_fbs, gv_total_dofs, total_dofs
    integer(kind=8) :: jmatt, npg, jcret
    aster_logical :: lMatr, lVect, lSigm, lVari, matsym
    character(len=4), parameter :: fami = 'RIGI'
    real(kind=8) :: rhs(MSIZE_TDOFS_MIX)
    type(HHO_matrix) :: lhs
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
! --- Get element parameters
!
    call elrefe_info(fami=fami, npg=npg)
!
! --- Number of dofs
    call hhoMecaDofs(hhoCell, hhoData, mk_cbs, mk_fbs, mk_total_dofs)
    call hhoTherDofs(hhoCell, hhoData, gv_cbs, gv_fbs, gv_total_dofs)
    total_dofs = mk_total_dofs+gv_total_dofs+gv_cbs
!
    if (nomopt /= "RIGI_MECA_TANG" .and. &
        nomopt /= "FULL_MECA" .and. &
        nomopt /= "FORC_NODA" .and. &
        nomopt /= "RAPH_MECA") then
        ASSERT(ASTER_FALSE)
    end if
!
! --- Properties of behaviour
!
    call hhoComporState%initialize(fami, nomopt, hhoCell%ndim, hhoCell%barycenter)
    hhoComporState%typmod(2) = 'GRADVARI'
!
! --- Initialize quadrature for the rigidity
!
    call hhoQuadCellRigi%initCell(hhoCell, npg)
!
! --- Initialize displacement, vari, ...
!
    call hhoMecaState%initialize(hhoCell, hhoData, hhoComporState)
    call hhoGVState%initialize(hhoCell, hhoData, hhoComporState)
!
! --- Compute Operators
!
    call hhoCalcOpGv(hhoCell, hhoData, hhoComporState%l_largestrain, &
                     hhoMecaState, hhoGvState)
!
! --- Compute local contribution
!
    call hhoGradVariLC(hhoCell, hhoData, hhoQuadCellRigi, &
                       hhoMecaState, hhoComporState, hhoGVState, lhs, rhs)
!
    call behaviourOption(nomopt, hhoComporState%compor, lMatr, lVect, lVari, lSigm)
!
! --- Save return code
!
    if (lSigm) then
        call jevech('PCODRET', 'E', jcret)
        zi(jcret) = hhoComporState%codret
    end if
!
! --- Save rhs
!
    if (lVect) then
        call writeVector('PVECTUR', total_dofs, rhs)
    end if
!
! --- Save of lhs
!
    if (lMatr) then
        call nmtstm(hhoComporState%carcri, jmatt, matsym)
!
        if (matsym) then
            call lhs%write('PMATUUR', ASTER_TRUE)
        else
            call lhs%write('PMATUNS', ASTER_FALSE)
        end if
    end if
!
    call lhs%free()
    call hhoMecaState%free()
    call hhoGVState%free()
!
end subroutine
