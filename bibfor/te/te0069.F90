! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine te0069(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
    use FE_eval_module
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/nlcomp.h"
    character(len=16) :: option, nomte
!
!     BUT:
!       CALCUL DES FLUX DE TEMPERATURE AUX POINTS DE GAUSS
!       OPTION : 'FLUX_ELGA'
!
! ---------------------------------------------------------------------
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer :: icamas, kp, imate, itemps
    integer :: icodre(1)
    character(len=32) :: phenom
    real(kind=8) :: time, Kglo(3, 3), fluglo(3), dtpg(3), tpg
    real(kind=8), pointer :: flux(:) => null()
    real(kind=8), pointer :: tempi(:) => null()
!
! ----------------------------------------------------------------------
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PTEMPER', 'L', vr=tempi)
    call jevech('PFLUXPG', 'E', vr=flux)
!
    time = zr(itemps)
!
    call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
    if (phenom .eq. 'THER_ORTH') then
        call jevech('PCAMASS', 'L', icamas)
    end if
!
    do kp = 1, FEQuadCell%nbQuadPoints
        tpg = FEEvalFuncScal(FEBasis, tempi, FEQuadCell%points_param(1:3, kp))
        dtpg = FEEvalGradVec(FEBasis, tempi, FEQuadCell%points_param(1:3, kp))
        call nlcomp(phenom, imate, icamas, FECell%ndim, FEQuadCell%points(1:3, kp), time, &
                    tpg, Kglo)
        fluglo = matmul(Kglo, dtpg)
        flux(FECell%ndim*(kp-1)+1:FECell%ndim*(kp-1)+FECell%ndim) = -fluglo(1:FECell%ndim)
    end do
!
end subroutine
