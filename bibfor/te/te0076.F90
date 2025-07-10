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
subroutine te0076(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
!
    implicit none
#include "jeveux.h"
#include "asterfort/nlcomp.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:   OPTION : 'RIGI_THER'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! ......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    character(len=8), parameter :: famiR = "RIGI"
    integer(kind=8) :: icamas, kp, imate, itemps
    integer(kind=8) :: icodre(1)
    character(len=32) :: phenom
    real(kind=8) ::   time
    real(kind=8) :: rigi(MAX_BS, MAX_BS)
    real(kind=8) ::  valQPK(3, 3, MAX_QP)
! ----------------------------------------------------------------------
    call FECell%init()
    call FEQuadCell%initCell(FECell, famiR)
    call FEBasis%initCell(FECell)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTR', 'L', itemps)
    time = zr(itemps)
!
    call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
!   pour stopper le calcul si PCAMASS n'est pas disponible
    if (phenom .eq. 'THER_ORTH') then
        call jevech('PCAMASS', 'L', icamas)
    end if
!
    valQPK = 0.d0
    do kp = 1, FEQuadCell%nbQuadPoints
        call nlcomp(phenom, famiR, kp, imate, FECell%ndim, FEQuadCell%points(1:3, kp), &
                    time, 0.d0, valQPK(1:3, 1:3, kp))
    end do
!
    call FEStiffJacoScal(FEQuadCell, FEBasis, valQPK, rigi)
    call writeMatrix("PMATTTR", FEBasis%size, FEBasis%size, ASTER_TRUE, rigi)
!
end subroutine
