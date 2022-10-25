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

subroutine te0449(nomopt, nomte)
!
use HHO_type
use HHO_utils_module
use HHO_size_module
use HHO_quadrature_module
use HHO_ther_module
use HHO_basis_module
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
#include "asterfort/rccoma.h"
#include "asterfort/utmess.h"
#include "asterfort/rcvalb.h"
#include "asterfort/lteatt.h"
#include "jeveux.h"
#include "blas/dsyr.h"

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
    type(HHO_basis_cell) :: hhoBasisCell
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    integer :: cbs, fbs, total_dofs, npg, ipg, imate, itemps
    integer :: icod(1)
    character(len=4) :: fami
    character(len=8) :: typmod(2)
    character(len=32) :: phenom
    real(kind=8) :: coeff, cp(1), time
    real(kind=8), dimension(MSIZE_CELL_SCAL):: BSCEval
    real(kind=8), dimension(MSIZE_TDOFS_SCAL, MSIZE_TDOFS_SCAL) :: lhs
!
! --- Get HHO informations
!
    call hhoInfoInitCell(hhoCell, hhoData)
!
! --- Get element parameters
!
    fami = 'MASS'
    call elrefe_info(fami=fami, npg=npg)
!
! --- Number of dofs
    call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    ASSERT(cbs <= MSIZE_CELL_SCAL)
    ASSERT(fbs <= MSIZE_FACE_SCAL)
    ASSERT(total_dofs <= MSIZE_TDOFS_SCAL)
!
    if (nomopt /= "MASS_THER") then
        ASSERT(ASTER_FALSE)
    end if
!
! --- Initialize quadrature for the mass
!
    call hhoQuadCellMass%initCell(hhoCell, npg)
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
    if (lteatt('LUMPE','OUI')) then
        ASSERT(ASTER_FALSE)
    endif
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPSR', 'L', itemps)
!
    call rccoma(zi(imate), 'THER', 1, phenom, icod(1))
    time = zr(itemps)
!
! --- Compute local contribution
!
! ----- init basis
    call hhoBasisCell%initialize(hhoCell)
    lhs = 0.d0
!
! ----- Loop on quadrature point
!
    do ipg = 1, hhoQuadCellMass%nbQuadPoints
! --------- Get rho * Cp
!
        if (phenom .eq. 'THER') then
            call rcvalb(fami, ipg, 1, '+', zi(imate),&
                        ' ', phenom, 1, 'INST', [time],&
                        1, 'RHO_CP', cp, icod(1), 1)
        else if (phenom .eq. 'THER_ORTH') then
            call rcvalb(fami, ipg, 1, '+', zi(imate),&
                        ' ', phenom, 1, 'INST', [time],&
                        1, 'RHO_CP', cp, icod(1), 1)
        else
            call utmess('F', 'ELEMENTS2_63')
        endif
!
! --------- Eval bais function at the quadrature point
!
        call hhoBasisCell%BSEval(hhoCell, hhoQuadCellMass%points(1:3,ipg), 0, &
                            hhoData%cell_degree(), BSCEval)
!
! --------  Eval massMat
!
        coeff = cp(1) * hhoQuadCellMass%weights(ipg)
        call dsyr('U', cbs, coeff, BSCEval, 1, lhs, MSIZE_TDOFS_SCAL)
    end do
!
! ----- Copy the lower part
!
    call hhoCopySymPartMat('U', lhs(1:cbs, 1:cbs))
!
! --- Save lhs
!
    call hhoRenumTherMat(hhoCell, hhoData, lhs)
    call writeMatrix('PMATTTR', total_dofs, total_dofs, ASTER_TRUE, lhs)
!
end subroutine
